import hashlib
import random


class LamportKeyGen:
    def __init__(self, seed):
        self.seed = seed
        self.step = -1
        self.reload_generators()
        self.X = []
        self.Y = []

    def reload_generators(self):
        self.step = -1
        random.seed(self.seed)

    def generate_private_key(self, fast=False):
        if fast:
            random.getrandbits(256)
            random.getrandbits(256)
            return

        self.X = []
        for i in range(256):
            x1 = hex(random.getrandbits(256))[2:].zfill(64)
            w1 = hex(random.getrandbits(256))[2:].zfill(64)
            self.X.append((x1, w1))

    def generate_public_key(self, fast=False):
        if fast:
            return

        self.Y = []
        for i in range(256):
            h_x1 = hashlib.sha256(self.X[i][0].encode('utf-8')).hexdigest()
            h_w1 = hashlib.sha256(self.X[i][1].encode('utf-8')).hexdigest()
            self.Y.append([h_x1, h_w1])

    def signature(self, msg):
        M_hex = hashlib.sha256(msg).hexdigest()
        M_bin = bin(int(M_hex, 16))[2:].zfill(256)
        sign = ''
        for i, bit in enumerate(M_bin):
            sign += self.X[i][int(bit)]
        return sign

    def get_public_hash(self):
        h = hashlib.sha256()
        for pair in self.Y:
            h.update(pair[0].encode('utf-8'))
            h.update(pair[1].encode('utf-8'))
        return h.hexdigest().zfill(64)

    def generate_new_keys(self, index=-1, fast=False):
        if index == -1:
            self.step += 1
            self.generate_private_key(fast=fast)
            self.generate_public_key(fast=fast)
        else:
            if index < self.step:
                self.reload_generators()
            if index == self.step:
                return
            else:
                for _ in range(index - self.step - 1):
                    self.generate_new_keys(fast=fast)
                self.generate_new_keys()