#include <iostream>
#include "GaloisField.h"
#include <cassert>
#include <thread>
#include <mutex>
#include "McEliceCryptoC.h"

constexpr const uint32_t GoppaDegree = 8;
int main() {


    auto* myCipher = static_cast<GaloisMath::LinearGoppaCodeGen<8, 13, 299> *>(createCipher());

    std::string pubKey;
    std::cout << "Key size: " << getPublicKeySize(myCipher) << std::endl;
    pubKey.resize(getPublicKeySize(myCipher));
    getPublicKey(myCipher, pubKey.data());
//    std::cout << "Pubkey: " << pubKey << std::endl;
    printf(pubKey.c_str());
    std::string message = "Hello!!";
    std::cout << "Initial message: " << message << std::endl;

    message.resize(getDecryptedMessageSize(myCipher));
    std::string encrypted;
    encrypted.resize(getEncryptedMessageSize(myCipher));

    encrypt(message.c_str(), pubKey.c_str(), encrypted.data());
    std::cout << "Encrypted message: " << encrypted << std::endl;
    std::string message_back;
    message_back.resize(getDecryptedMessageSize(myCipher));
    decrypt(myCipher, encrypted.c_str(), message_back.data());
    std::cout << "Decrypted message: " << message_back << std::endl;

    freeCipher(myCipher);

    return 0;
}
