f = open ('centroids', 'rt'); 

while 1==1 :
	x1, x2, x3 = [float (x) for x in f.readline().split()]; 
	y1, y2, y3 = [float (x) for x in f.readline().split()];
	d = ((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2)+(x3-y3)*(x3-y3))**0.5; 
	print '{0:2f}'.format(d);
