import os, socket, time

BUFFERSIZE = 512

class pTemp:
    def __init__(self,location, name, logFile, seed, basePort):
        os.system('touch '+logFile+'.txt')
        portID = basePort
        while (portID<basePort+20*10):
            try:
                cmd = 'xterm -e "./'+location+'/process '+name+' '+logFile+'.txt '+str(seed)+' '+str(portID)+'" &'
                os.system(cmd)
                self.s = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
                self.s.bind(('127.0.0.1', portID))
                #self.s.bind(('143.167.50.72', portID))
                break
            except:
                portID+=20
        print cmd
        self.s.listen(100)
        self.conn, self.addr = self.s.accept()
        print 'Connection:', self.addr

    def stream(self,c):
		self.data = self.conn.recv(BUFFERSIZE)
		if self.data:
			self.conn.send(str(c))
		#time.sleep(0.005)

    def out(self):
        return self.data.split(',')[:-1]

    def quit(self):
        self.stream(0)
        self.conn.close()

    def step(self, n=1):
        for t in range(n):
            self.stream(1)
        return self.data

    def show(self):
        self.stream(2)

    def hide(self):
        self.stream(3)

    def setv(self,stateI,stateJ,val):
        self.stream('4,'+str(stateI)+','+str(stateJ)+','+str(val))

    def save(self, fileName):
        self.stream('5,'+fileName)

    def load(self, fileName):
        self.stream('6,'+fileName)

    def join(self, portID):
        self.stream('7,'+str(portID))

    def misc(self):
        self.stream('8')


class pBase:
    def __init__(self, portID):
        self.s = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        self.s.bind(('127.0.0.1', portID))
        self.s.listen(1)
        self.conn, self.addr = self.s.accept()
        print 'Connection:', self.addr

    def stream(self,c):
		self.data = self.conn.recv(BUFFERSIZE)#note has been 10000
		if self.data:
			self.conn.send(str(c))
		#time.sleep(0.005)

    def quit(self):
        self.conn.close()


'''
    BESPOKE WRAPPER CLASSES FOR EACH PROCESS
'''
