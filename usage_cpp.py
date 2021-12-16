from ctypes import *

lib = cdll.LoadLibrary('./McEliece/libMcElieceCryptoLib.so')
createCipher = lib.createCipher
getPublicKeySize = lib.getPublicKeySize  # (gen)
getPublicKey = lib.getPublicKey  # (gen, pubKey)
# getPublicKey.argtypes = c_void_p, c_char_p
getEncryptedMessageSize = lib.getEncryptedMessageSize
getDecryptedMessageSize = lib.getDecryptedMessageSize


_encrypt = lib.encrypt
_decrypt = lib.decrypt


class McEliece:
    def __init__(self):
        self.gen = createCipher()
        self.keySize = getPublicKeySize()
        self.pubKey = create_string_buffer(self.keySize)
        getPublicKey(self.gen, self.pubKey)
        self.enc_size = getEncryptedMessageSize(self.gen)
        self.dec_size = getDecryptedMessageSize(self.gen)
        self.decrypted = create_string_buffer(self.dec_size)
        self.encrypted = create_string_buffer(self.enc_size)
        self.alien_encrypted = create_string_buffer(self.enc_size)

    def set_pubKey(self, alien_pubKey):
        self.pubKey.value = alien_pubKey

    def encrypt(self, msg):
        _encrypt(msg.encode('utf-8'), self.pubKey, self.encrypted)
        return self.encrypted.value

    def decrypt(self, alien_msg):
        self.alien_encrypted.value = alien_msg
        _decrypt(self.gen, self.alien_encrypted, self.decrypted)
        return self.decrypted.value

if __name__ == '__main__':
    mc = McEliece()
    key = mc.pubKey.raw
    mc.pubKey.raw = key
    print(key)
    a = mc.encrypt('hello')
    print(a)
    print(mc.decrypt(a))
