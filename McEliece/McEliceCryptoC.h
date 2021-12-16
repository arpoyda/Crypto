//
// Created by Бушев Дмитрий on 14.12.2021.
//

#ifndef MCELIECECRYPTO_MCELICECRYPTOC_H
#define MCELIECECRYPTO_MCELICECRYPTOC_H

extern "C" void* createCipher();

extern "C" void freeCipher(void* cipher);

extern "C" uint32_t getPublicKeySize(void* codeGen);

extern "C" uint32_t getDecryptedMessageSize(void* codeGen);
extern "C" uint32_t getEncryptedMessageSize(void* codeGen);

extern "C" void getPublicKey(void* codeGen, char* key);

extern "C" void decrypt(void* codeGen, const char* encryptedMessage, char* decryptedMessage);

extern "C" void encrypt(const char* message, const char* key, char* encryptedMessage);

#endif //MCELIECECRYPTO_MCELICECRYPTOC_H
