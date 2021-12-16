import socket
import secrets
from MerkleSign import MerkleSign
from MerkleDesign import MerkleDesign
from usage_cpp import McEliece


class ClientError(Exception):
    pass


class ServerError(Exception):
    pass


class Client:
    def __init__(self, host='127.0.0.1', port=62000, timeout=None, buffer_size=2**16, message_number=8, log=False):
        self.host = host
        self.port = port
        self.timeout = timeout
        self.mceliece = McEliece()
        self.buffer_size = buffer_size
        print('[Client]: Waiting for Merkle establishing...')
        self.merkle_sign = MerkleSign(message_number, secrets.randbits(256))
        print('[Client]: Merkle established!')

        try:
            self.connection = socket.create_connection((host, port), timeout)
        except socket.error as err:
            raise ClientError("[Client]: Cannot create connection", err)

    def _send(self, data):
        try:
            self.connection.sendall(data)
        except socket.error as err:
            raise ClientError("Error sending data to server", err)

    def _read(self):
        try:
            data = self.connection.recv(self.buffer_size)
        except socket.error as err:
            raise ClientError("Error reading data from socket", err)
        return data

    def send_root(self):
        self._send(self.merkle_sign.root.encode())
        response = self._read()
        self.mceliece.pubKey.raw = response
        print('[Client]: McEliece public key received!')

    def send_msg(self, msg):
        responses = []
        b, Y, Auth = self.merkle_sign.sign_message(msg.encode('utf-8'))

        self._send(self.mceliece.encrypt(msg))
        responses.append(self._read().decode() == 'OK')

        self._send(b.encode('utf-8'))
        responses.append(self._read().decode() == 'OK')

        self._send(Y.encode('utf-8'))
        responses.append(self._read().decode() == 'OK')

        self._send(Auth.encode('utf-8'))
        responses.append(self._read().decode() == 'OK')

        if not all(responses):
            print('Some problems with responses from Server')

    def run(self):
        self.send_root()
        while True:
            msg = input() + '\n'
            self.send_msg(msg)
            if msg.lower() == 'bye':
                break


class Server:
    def __init__(self, addr='127.0.0.1', port=62000, buffer_size=2**16, log=False):
        self.merkle_design = None
        self.mceliece = McEliece()
        self.buffer_size = buffer_size
        self.socket_server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket_server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            self.socket_server.bind((addr, port))
            self.socket_server.listen()
            if log:
                print('[SERVER]: Socket server started and listening at port {}:'.format(port))
        except socket.error as err:
            raise ServerError("Cannot create connection", err)

    def run(self):
        while True:
            conn, addr = self.socket_server.accept()
            with conn:
                print(f'[SERVER]: Connected: {addr}')
                self.merkle_design = MerkleDesign(conn.recv(self.buffer_size).decode('utf-8'))
                conn.sendall(self.mceliece.pubKey.raw)
                print('[SERVER]: Root received!')
                while True:
                    msg_encrypted = conn.recv(self.buffer_size)
                    msg = self.mceliece.decrypt(msg_encrypted).decode()
                    print('McEliece decrypted correctly!')
                    if not msg:
                        break
                    conn.sendall('OK'.encode())
                    b = conn.recv(self.buffer_size).decode()
                    conn.sendall('OK'.encode())
                    Y = conn.recv(self.buffer_size).decode()
                    conn.sendall('OK'.encode())
                    Auth = conn.recv(self.buffer_size).decode()
                    conn.sendall('OK'.encode())

                    if self.merkle_design(msg.encode('utf-8'), b, Y, Auth):
                        print('[SERVER]: Message:', msg)
                        if msg.lower() == 'bye':
                            break
                break




if __name__ == '__main__':
    pass
