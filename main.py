from socket_connection import Client, Server
from threading import Thread
import time

my_addr = input('Enter your IP and PORT: ')
friend_addr = input('Enter friend\'s IP and PORT: ')

my_addr = my_addr if my_addr else '127.0.0.1 62000'
friend_addr = friend_addr if friend_addr else '127.0.0.1 62000'

# my_addr = my_addr if my_addr else '25.89.157.192 62000'
# friend_addr = friend_addr if friend_addr else '25.25.183.48 62000'

s = Server(my_addr.split()[0], int(my_addr.split()[1]))
t_s = Thread(target=s.run)
t_s.start()

input('>>> Print something to create CONNECTION: ')

c = Client(friend_addr.split()[0], int(friend_addr.split()[1]))
t_c = Thread(target=c.run)
t_c.start()

t_s.join()

#Enter your IP and PORT: 25.89.157.192 62000
#Enter friend's IP and PORT: 25.25.183.48 62000