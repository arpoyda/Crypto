import hashlib
from textwrap import wrap
from MerkleSign import xor_hex


class MerkleDesign:
    def __init__(self, public_key):
        self.public_key = public_key

    def __call__(self, msg, b, Y, Auth):
        Y_buff = wrap(Y, 64)
        Y1 = Y_buff[::2]
        Y2 = Y_buff[1::2]
        Y = list(zip(Y1, Y2))

        Auth = wrap(Auth, 64)

        return self.lamport_signature_check(msg, b, Y) & self.merkle_signature_check(Y, Auth)

    def merkle_signature_check(self, Y, Auth):
        h = hashlib.sha256()
        for pair in Y:
            h.update(pair[0].encode('utf-8'))
            h.update(pair[1].encode('utf-8'))
        H_ij = h.hexdigest()
        for H_auth in reversed(Auth):
            H_ij = xor_hex(H_ij, H_auth)

        if H_ij == self.public_key:
            print('Merkle signature is MATCHED!')
            return True
        else:
            print('Merkle signature is MISMATCHED!')
            return False

    def lamport_signature_check(self, msg, b, Y):
        hashed_msg = self.hash_message(msg)
        hashed_b = self.hash_b(b)
        hashed_Y = self.concat_Y(Y, hashed_msg)
        if hashed_b == hashed_Y:
            print('Lamport signature is MATCHED!')
            return True
        else:
            print('Lamport signature is MISMATCHED!')
            return False

    def hash_message(self, msg):
        M_hex = hashlib.sha256(msg).hexdigest()
        M_bin = bin(int(M_hex, 16))[2:].zfill(256)
        return M_bin

    def hash_b(self, b):
        hex_length = int(256 / 4)
        b_hashed = ''
        for i in range(256):
            b_hashed += hashlib.sha256(b[i * hex_length: (i + 1) * hex_length].encode('utf-8')).hexdigest()
        return b_hashed

    def concat_Y(self, Y, hash_msg):
        concated = ''
        for i, bit in enumerate(hash_msg):
            concated += Y[i][int(bit)]
        return concated