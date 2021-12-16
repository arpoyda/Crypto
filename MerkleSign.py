import hashlib

from Lamport import LamportKeyGen
from MerkleCache import Cache

def xor_hex(x1, x2):
    return hex(int(x1, 16) ^ int(x2, 16))[2:].zfill(64)

class MerkleSign:
    def __init__(self, N, seed):
        self.N = N  # 2^N messages can be signed
        self.keygen = LamportKeyGen(seed)
        self.root = self.get_public_root(0)
        self.cache = Cache(N)
        self.step = 0
        self.TEST_counter = 0

    def get_public_root(self, layer):  # TODO: change to recalc_or_get because it is equal
        if layer == 0:
            self.keygen.reload_generators()

        if layer == self.N:
            self.keygen.generate_new_keys()
            return self.keygen.get_public_hash()

        return xor_hex(self.get_public_root(layer + 1), self.get_public_root(layer + 1))

    def recalc_or_get(self, index, layer):
        if self.cache[index, layer]:
            return self.cache[index, layer]

        if layer == self.N:
            self.keygen.generate_new_keys(index)
            return self.keygen.get_public_hash()

        return xor_hex(self.recalc_or_get((index << 1) + 0, layer + 1), self.recalc_or_get((index << 1) + 1, layer + 1))

    def find_auth_path(self):
        auth_path = []
        for i in range(self.N):
            shifted_index = self.step >> i
            inv_bit = (shifted_index & 1) ^ 1
            auth_path.append((shifted_index >> 1 << 1) + inv_bit)

        return list(reversed(auth_path))

    def fill_cache(self):
        auth_path = self.find_auth_path()
        for i, index in enumerate(auth_path):
            layer = i + 1
            self.cache[index, layer] = self.recalc_or_get(index, layer)  # excess copy BUT WE ARE IN PYTHON

    def sign_message(self, msg):
        self.fill_cache()
        self.keygen.generate_new_keys(index=self.step)
        b = self.keygen.signature(msg)
        Y = ''.join(sum(self.keygen.Y, []))
        Auth = ''.join(self.cache.get_auth_path())
        self.step += 1
        return b, Y, Auth

