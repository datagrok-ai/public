import psutil
from flask import Flask
import json
app = Flask(__name__)

 # CPU load info:
def cpu ():
	cpu_per = int (psutil.cpu_percent ()) 
	return cpu_per
 # Memory load info:
def mem ():
 # mem = psutil.virtual_memory () #memory full information;
	mem_total = int(psutil.virtual_memory()[0]/1024/1024)
	mem_used = int(psutil.virtual_memory()[3] / 1024 / 1024)
	mem_per = int(psutil.virtual_memory()[2])
	mem_info = {
		'mem_total' : mem_total,
		'mem_used' : mem_used,
		'mem_per' : mem_per
	}
	return mem_info
 # HDD/SSD info:
def disk ():
	d_per = int(psutil.disk_usage('/')[3]) 
	
	disk_info = {
		'used' : d_per,
		
	}
	return disk_info
 # Network traffic monitor:
def network ():
 # network = psutil.net_io_counters () # all info;
	network_sent = int (psutil.net_io_counters () [0] / 8/1024) #kb per second
	network_recv = int(psutil.net_io_counters()[1]/8/1024)
	network_info = {
		'network_sent' : network_sent,
		'network_recv' : network_recv
	}
	return network_info

	
@app.route('/')
def root():
	return json.dumps({
		'CPU_average_load, % ': cpu(), 
		'Memory_usage, MB ': mem(),
		'Disk_usage, %' : disk(),
		'Network_usage, kb/sec' : network()
	})

#main()
if __name__ == '__main__':
    app.run(host="0.0.0.0")
