//
// Created by Бушев Дмитрий on 04.12.2021.
//

#include "GaloisField.h"
#include <stack>
#include <iostream>
#include "McEliceCryptoC.h"


namespace GaloisMath{

    RandomEngine randomEngine;

    uint32_t pow(uint32_t x, uint32_t degree){
        if(degree == 0)
            return 1u;
        if(degree == 1)
            return x;
        uint32_t lhs = degree / 2u;
        uint32_t rhs = degree - lhs;

        return pow(x, lhs) * pow(x, rhs);
    }




    std::pair<uint32_t, uint32_t> divide_with_remainder(uint32_t num, uint32_t divisor){
        std::pair<uint32_t, uint32_t> ret;
        ret.first = num / divisor;
        ret.second = num - divisor * ret.first;
        return ret;
    }

    uint32_t inverse_modulo(uint32_t a, uint32_t modulo){
        if(a == 1u)
            return 1u;

        std::vector<std::pair<uint32_t, uint32_t>> euclid_stack;
        a = a % modulo;

        // Direct euclid algorithm
        uint32_t greater = modulo, smaller = a;
        euclid_stack.emplace_back(0u, greater);
        euclid_stack.emplace_back(0u, smaller);
        while(smaller != 1){
            uint32_t div = greater / smaller;
            uint32_t rem = greater % smaller;
            euclid_stack.emplace_back(div, rem);

            greater = smaller;
            smaller = rem;
        }

        // Reverse euclid algorithm

        // q * ac + p * bc = 1
        int q, ac, p, bc;

        q = greater;
        p = smaller;
        ac = 1;
        bc = -euclid_stack.back().first;



        for(int i = euclid_stack.size() - 1; i >= 3; --i){

            int tmp_ac = ac;
            ac = bc;
            bc = tmp_ac - euclid_stack.at(i - 1).first * bc;
            p = q;
            q = euclid_stack.at(i - 3).second;

        }

        if(bc < 0)
            bc += modulo;

        bc %= modulo;


        return bc;
    }

    uint32_t divide_modulo(uint32_t a, uint32_t b, uint32_t modulo){
        return (a * inverse_modulo(b, modulo)) % modulo;
    }


}
using GeneralCodeGen = GaloisMath::LinearGoppaCodeGen<8, 13, 0x12B>;

extern "C" void* createCipher(){
    return new GeneralCodeGen{GeneralCodeGen::GoppaCoefs{53,100,17,229,248,45,120,152,113,131,133,197,103,129}};
}

extern "C" void freeCipher(void* cipher){
    delete reinterpret_cast<GeneralCodeGen*>(cipher);
}

extern "C" uint32_t getPublicKeySize(void* codeGen){
    return 3 * sizeof(uint32_t) + 256 * 19 * 8 / 8;
}

extern "C" uint32_t getDecryptedMessageSize(void* codeGen){
    return 19;
}
extern "C" uint32_t getEncryptedMessageSize(void* codeGen){
    return 32;
}

extern "C" void getPublicKey(void* codeGen, char* key){
    auto pubKey = reinterpret_cast<GeneralCodeGen*>(codeGen)->getPublicKey();
    for(uint32_t i = 0; i < 256; ++i)
        for(uint32_t j = 0; j < 19 * 8; ++j) {
            if (pubKey.first.at(i, j) != 0u)
                key[(j * 256 + i) / 8 + 3 * sizeof(uint32_t)] =
                        key[(j * 256 + i) / 8 + 3 * sizeof(uint32_t)] | (1u << ((j * 256 + i) % 8));
        }

    *(reinterpret_cast<uint32_t *>(key)) = 256;
    *(reinterpret_cast<uint32_t *>(key + sizeof(uint32_t))) = 19 * 8;
    *(reinterpret_cast<uint32_t *>(key + 2 * sizeof(uint32_t))) = 13;


}

extern "C" void decrypt(void* codeGen, const char* encryptedMessage, char* decryptedMessage){
    GaloisMath::Matrix<256, 1, GaloisMath::GaloisElem<2u>> encMsg;
    for(uint32_t i = 0; i < 32; ++i){
        for(uint32_t j = 0; j < 8; ++j){
            encMsg.at(i * 8 + j, 0) = (encryptedMessage[i] & (1u << j)) != 0 ? 1u : 0u;
        }
    }

    auto decrypted = reinterpret_cast<GeneralCodeGen*>(codeGen)->decrypt(encMsg);
    for(uint32_t i = 0; i < 19; ++i){
        decryptedMessage[i] = '\0';
        for(uint32_t j = 0; j < 8; ++j){
            if(decrypted.at(i * 8 + j, 0) != 0)
                decryptedMessage[i] = decryptedMessage[i] | (1u << j);
        }
    }

}
extern "C" void encrypt(const char* message, const char* key, char* encryptedMessage){
    GaloisMath::Matrix<19 * 8, 1, GaloisMath::GaloisElem<2u>> msg;

    std::pair<GaloisMath::Matrix<256, 19 * 8, GaloisMath::GaloisElem<2u>>, uint32_t> key_;

    if(*reinterpret_cast<const uint32_t*>(key) != 256 ||
       *reinterpret_cast<const uint32_t*>(key + sizeof(uint32_t)) != 19 * 8)
        throw std::runtime_error("Broken key");

    key_.second = *reinterpret_cast<const uint32_t*>(key + 2 * sizeof(uint32_t));

    for(uint32_t i = 0; i < 256; ++i)
        for(uint32_t j = 0; j < 19 * 8; ++j) {
            key_.first.at(i, j) = (key[(j * 256 + i) / 8 + 3 * sizeof(uint32_t)] & (1u << ((j * 256 + i) % 8))) != 0u ? 1u : 0u;
        }
    for(uint32_t i = 0; i < 19 * 8; ++i){
        msg.at(i, 0) = (message[i / 8] & (1u << (i % 8))) != 0u ? 1u : 0u;
    }


    auto encrypt = GaloisMath::McElieceCipher<256, 19 * 8>{key_}.encrypt(msg);

    for(uint32_t i = 0; i < 32; ++i) {
        encryptedMessage[i] = '\0';
        for (uint32_t j = 0; j < 8; ++j) {
            if (encrypt.at(i * 8 + j, 0) != 0u)
                encryptedMessage[i] = encryptedMessage[i] | (1u << j);
        }
    }
}