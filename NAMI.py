# Copyright (c) 2014, Durham University
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the Durham University nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL DURHAM UNIVERSITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Please cite: Groftehauge et al. 2014. Acta Cryst D.

from __future__ import unicode_literals
import sys, os, numpy, csv, math, re, ntpath, itertools, decimal 
from functools import partial
from collections import Counter

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import pyqtSlot, SIGNAL, SLOT
from PyQt4.Qt import *
#from PyQt4.Qt import WindowStaysOnTopHint

#from numpy import genfromtxt, NaN, Inf, arange, isscalar, array


# MATPLOTLIB

#need to check whether I neeed this one, don't dare to remove it yet
#import matplotlib as mpl

#from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
# These are needed to plot the waterfall plot, which is a hack of the polygon plots. 
# See: http://matplotlib.org/examples/mplot3d/polys3d_demo.html
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.backend_bases import FigureManagerBase, KeyEvent, MouseEvent
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import colorConverter
from matplotlib.font_manager import FontProperties
import matplotlib.colors as col
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

# SCIPY
from scipy.interpolate import interp1d
import scipy.signal
from pylab import * #imports matplotlib as mpl and numpy as np and makes it behave like MatLab
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy import stats 


# PRETTYPLOTLIB
#import prettyplotlib as ax
#from prettyplotlib import brewer2mpl

###################################################################################################################################

progname = os.path.basename(sys.argv[0])
progversion = "0.1"

class Dictionaries():

	def namedictionary(self):
		''' creates a dictionary that relates well  number to well names'''

		# WELL NAMES 
		labelc = ['1', '2', '3', '4', '5', '6', '7', '8', '9',
		          '10', '11', '12']
		labelr = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H' ]
		names = []
		for c in labelc:
			for r in labelr:
				names.append(r+c)
		convert = lambda text: int(text) if text.isdigit() else text
		alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
		a = sorted(names, key = alphanum_key)

		# WELL NUMBERS 
		wells = []
		for w in range(1,97):
			wells.append(('Well %s' % w))

		self.dict_numb_name = dict(zip(wells, a))
		return self.dict_numb_name


	def solutionsfile(self, *args):

		''' Reads in the solution file into a dictionary. The key name is the well number and the values are all the different variables in the solution file. '''

		try:
			args = list(args)[0]
			fsolution = numpy.genfromtxt(args, delimiter=',', dtype=None)

			adict = {}
			for w in range(1,97):
				if w == fsolution[w-1][0]:
					adict[('Well %s' % (w))] = (round(fsolution[w-1][1],2),round(fsolution[w-1][2],4),fsolution[w-1][3],fsolution[w-1][4],fsolution[w-1][5],round(fsolution[w-1][6],2),int(fsolution[w-1][7]),int(fsolution[w-1][8]))
			self.dict_sol_numb = adict
			#print adict['Well 44']	

		except:
			self.dict_sol_numb = {}
		return self.dict_sol_numb



class ApplicationWindow(QtGui.QMainWindow):
	''' This is the main body of the program: contains functions for doing the calculations as well as the main window of the GUI. Coordinates signal transfer between classes (and other windows)'''

	def __init__(self):
		QtGui.QMainWindow.__init__(self)
		self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
		#self.resize()
		self.create_main_frame()
		self._generator = None
		self._timerId = None
		self.th_data = []

	def create_main_frame(self):
		
		#self.main_frame = QWidget()
		self.main_widget = QtGui.QWidget(self)
		#self.setFixedSize(950, 800)



		##############################################

		# FOR THE TABLE 
		self.widgetT = TablePopup()
		self.widgetT.setGeometry(QRect(100, 100, 400, 300))

		self.connect(self.widgetT, SIGNAL('mySigT'), self.tableParameters)
		self.connect(self.widgetT, SIGNAL('Table'), self.tempTable)


		#############################################

		# QUIT BAR
		self.file_menu = QtGui.QMenu('&File', self)
		self.file_menu.addAction('&Quit', self.fileQuit,
		                         QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
		self.menuBar().addMenu(self.file_menu)

		##############################################


		# OPEN BAR, in the old version they were only on the menu button
		self.open_menu = QtGui.QMenu('&Open', self)
		self.file_menu.addAction('&Open Dataset', self.fileOpen)
		self.menuBar().addMenu(self.file_menu)


		self.open_menu = QtGui.QMenu('&Open', self)
		self.file_menu.addAction('&Open solutionsfile', self.openSolFile)
		self.menuBar().addMenu(self.file_menu)


		self.open_menu = QtGui.QMenu('&Open', self)
		self.file_menu.addAction('&Open Results', self.openResults)
		self.menuBar().addMenu(self.file_menu)



		##############################################


		# CALLING ANOTHER CLASS: THE PLOT MENU
		# this forms the window where we get the extra plots. it calls the class that contains the data handling
		self.plot_menu = QtGui.QMenu('&Data', self)
		self.menuBar().addSeparator()
		self.plot_menu.addAction('Plots for Analysis', self.plotMenu)
		self.plot_menu.addAction('Load table (from previous results)', self.widgetT.showTable)
		self.plot_menu.addAction('Re-input parameters', self.resultsinputwindow)
		self.menuBar().addMenu(self.plot_menu)
		self.widgetplot = AnalysisPlotPopup()
		self.widgetplot.setGeometry(400,200,950,900)


		############################################


		# HELP BAR
		# adding details about what to cite and how to use the program should be inserted here. (connect it to self.help_menu)
		self.help_menu = QtGui.QMenu('&Help', self)
		self.menuBar().addSeparator()
		self.menuBar().addMenu(self.help_menu)
			# ABOUT 
		self.help_menu.addAction('&About', self.about)



		
		##############################################

		# might want to make this into a gridlayout to be able to position the buttons more specifically 
		l = QtGui.QVBoxLayout(self.main_widget)


		#############################################


		#DATA ANALYSIS: buttons that coordinate into the input file parameters and the data analysis process.  

		# START LOOP 
		self.button1 = QtGui.QPushButton('Set Input Parameters')
		self.button1.clicked.connect(self.start)
		l.addWidget(self.button1)
		self.widget = InputParametersPopup()
		#self.widget.setGeometry(QRect(1000, 1000, 1000, 1000))
		# read in data emitted from the other class
		self.connect(self.widget, SIGNAL('mySig'), self.parameters)

		# CONTINUE LOOP 
		self.button2 = QtGui.QPushButton('Begin Data Analysis')
		self.button2.clicked.connect(self.continues)
		l.addWidget(self.button2)


		##############################################


		# FIGURE: Creating a canvas for matplotlib etc. 
		self.figure = plt.figure(dpi=100)
		self.canvas =  mpl.backends.backend_qt4agg.FigureCanvasQTAgg(self.figure)
		self.canvas.setMinimumSize(800,500)
		l.addWidget(self.canvas)
		self.canvas.setParent(self.main_widget) # main_widget parent instead of main_frame
		self.canvas.setFocusPolicy(Qt.ClickFocus)
		self.canvas.setFocus()


		########################################


		# Calling classes containing dictionaries of well names and numbers
		D = Dictionaries()
		self.b = D.namedictionary()


		########################################


		# NAVIGATION TOOLBAR
		# this is the little tiny toolbar at the bottom, where you can save individual figures etc. Don't think most buttons work on it though...
		self.mpl_toolbar = NavigationToolbar(self.canvas, self)
		l.addWidget(self.mpl_toolbar)
		self.canvas.mpl_connect('button_press_event', self.on_key_press)
		#self.main_widget.setFocus()
		self.setCentralWidget(self.main_widget)


		########################################

		# EXTRAS: Checkboxes, etc. 

		# CHECKBOX FOR AUTOMATIC MODE
		self.chick = QtGui.QCheckBox('Automatic', self)
		l.addWidget(self.chick)
		#self.connect(self.chick, QtCore.SIGNAL('toggled(bool)'), self.continues)

		# CHECKBOX FOR SAVEFILE 
		self.savfil = QtGui.QCheckBox('Save Graphs', self)
		l.addWidget(self.savfil)
		self.connect(self.savfil, QtCore.SIGNAL('toggled(bool)'), self.updategraph)

		# COMMENTS TO TXT FILE
		self.edit1 = QtGui.QLineEdit()
		l.addWidget(self.edit1)
		self.edit1.setPlaceholderText("Insert comments for the text file ")

		# STATUSBAR 
		self.statusBar().showMessage("All hail Queen Nelly!", 2000)

		###########################################


		# FOR SIGNALLING: CONNECTING BUTTONS
		self.connect(self, SIGNAL('datadict'), self.widgetplot.graphdata)
		self.connect(self, SIGNAL('colordict'), self.widgetplot.colouring)
		self.connect(self.widgetT, SIGNAL('mySigT'), self.widgetplot.values)
		self.connect(self, SIGNAL('colordict1'), self.widgetplot.colouring)
		self.connect(self, SIGNAL('results'), self.widgetplot.values)
		self.connect(self, SIGNAL('solvs'), self.widgetplot.canvastime)
		self.connect(self.widget, SIGNAL('opensolution'), self.openSolFile)
		self.connect(self.widget, SIGNAL('openrawdata'), self.fileOpen)
		self.connect(self.widget, SIGNAL('openresults'), self.openResults)
		self.connect(self, SIGNAL('opens'), self.widget.reciever_s)
		self.connect(self, SIGNAL('openr'), self.widget.reciever_r)
		self.connect(self, SIGNAL('open_rawdata'), self.widget.reciever_rawdata)

		########################################## 

	def fileQuit(self):
		# closes the window 
		self.close()

	def plotMenu(self):
		self.emit(QtCore.SIGNAL('solvs'), self.solvfilename)
		self.widgetplot.show()


	def on_key_press(self, event):

		''' This is for the table. It divides up the matplotwindow into squares where the contents of the solution file will show when you hover over it. For instance, where you to click and hover 
		over a normal denaturation plot, you would see the contents of a certain well. This is a slight bug, but I think its OK as long as no one notices. As you can see, it is a function that reacts to an event, it does not need to be passed any variables. The event that I have defined it to react to is obviously buttonpress.'''

		
		l_h = self.figure.get_window_extent().height
		l_w = self.figure.get_window_extent().width
		
		# we need to use those above to scale it 

		i_w = 800
		i_h = 426

		scale_w = i_w/l_w
		scale_h = i_h/l_h

		# need to fix these in case the table moves in pixels  
		width = 51/scale_w
		height = 35/scale_h
		start_x = 61/scale_w
		start_y = 376/scale_h
		coords = []

		# Building up the coordinate system, relating coordinates with the well number 
		for s in range(1,13,1):
			for d in range(1,9,1):
				coords.append((s,d)) 
		coords = [(x,y) for y,x in coords] # (r,c) structure
		grid = sorted(coords)
		adict = {}
		for w in range(1,97):
			adict[grid[w-1]] = ('Well %s' % (w))

		try:
			D = Dictionaries()
			sdict = D.solutionsfile(self.solvfilename)

		# populate an array 
			n = []

			for i in range(0,96):
				nh = grid[i][0]
				nw = grid[i][1]
				n.append([(start_x+(nw*width), start_y-(nh*height)), (start_x+((nw+1)*width), start_y-((nh+1)*height))])
				#if (start_x+(nw*width), start_y-(nh*height)) <= (event.x, event.y) <= ((start_x+((nw+1)*width), start_y-((nh+1)*height))):
				#if (start_x+(nw*width)) <= event.x <= (start_x+((nw+1)*width)):
				#print (start_y-(nh*height)), (start_y-((nh+1)*height))

			# detecting where the cursor is when the button is pressed 
			for i in range(len(n)):
				if n[i][0][0] < event.x < n[i][1][0]:
					if n[i][1][1] <= event.y <= n[i][0][1]:
						string = str(sdict[adict[grid[i]]]) # just showing contents 
						self.mpl_toolbar.set_message(str(adict[grid[i]]))
						self.canvas.setToolTip(string)
						print adict[grid[i]]

		except: # when no sol file has been submitted, it wont crash
			pass

		#button_press_handler(event, self.canvas, self.mpl_toolbar)
		
	def openSolFile(self):
		''' for opening the solutionfile'''

		self.solvfilename = str(QtGui.QFileDialog.getOpenFileName(self, 'Open File', os.curdir))
		self.emit(QtCore.SIGNAL('opens'), self.solvfilename) # emit solutionfilename so it can be shown in the window (when you have selected it)
	

		if self.solvfilename.endswith('.sol'): #making sure we only open .sol files
			pass
		else:
			QtGui.QMessageBox.critical(self, "Error", "It's not a .sol file, double check please.")


		return self.solvfilename

	def openResults(self):

		''' opening the resultsfile, similar to above function  '''

		self.resultsfile = str(QtGui.QFileDialog.getOpenFileName(self, 'Open File', os.curdir))
		self.emit(QtCore.SIGNAL('openr'), self.resultsfile)



		if self.resultsfile.endswith('.csv'):
			self.resultsinputwindow()
		else:
			QtGui.QMessageBox.critical(self, "Error", "That does not look like a .csv file.")

		print self.resultspopulation()
		self.emit(QtCore.SIGNAL('colordict1'), self.th_data)
		#self.emit(QtCore.SIGNAL('results'), float(self.referencet),float(self.sigfigchange))

	def resultsinputwindow(self): 		
		'''This is the window that asks the user for two values: the signficant change in temperature and the reference temperature. It is
		used when the user re-reads in old results and want to tweak the table/plots without having to restart the program. '''

		(text,truth)=QInputDialog.getText(self,"Get text","Reference temperature",QLineEdit.Normal,"50")
		(text1,truth1)=QInputDialog.getText(self,"Get text","Signficant change in temperature",QLineEdit.Normal,"2")
		if truth and truth1:
			self.referencet = text
			self.sigfigchange = text1
		else:
			pass

		# it then emits the signal afterwards 
		self.emit(QtCore.SIGNAL('results'), float(self.referencet),float(self.sigfigchange))

	def resultspopulation(self): # NEED TO MAKE THIS TAKE INTO CONSIDERATION IF USER SAID SOME OTHER VALUE 

		''' pretty sure this function processes the resultsfile into an array. '''

		if len(self.resultsfile) > 0:
			fresults = numpy.genfromtxt(self.resultsfile, delimiter=',', dtype=None)
			for w in range(3,99):
				if len(fresults[w][2]) == 0:
					self.th_data.append(float(fresults[w][1]))
				else:
					self.th_data.append((float(fresults[w][1]),float(fresults[w][2])))

			self.parsefile(int(fresults[1][0]),int(fresults[1][1]),int(fresults[1][2]), int(fresults[1][3]))
			#print "self, thdaa", self.th_data


	def fileOpen(self):
		self.openFile()

	def openFile(self):
		self.filename = str(QtGui.QFileDialog.getOpenFileName(self, 'Open File', os.curdir))
		self.emit(QtCore.SIGNAL('open_rawdata'), self.filename)


		if self.filename.endswith('.csv'):
			pass
		else:
			QtGui.QMessageBox.critical(self, "Error", "That does not look like a .csv file")
		return self.filename 

	def parameters(self, *args):
		''' Really stupid way of accepting parameters send from another class, in the end i never figured out the alternative '''
		self.args = list(args) # now its a list, so I can access the values easier
		return self.args

	def parsefile(self, column_temp=0, offset_temp=23, yval=3, columnwell=1, wells=97, increment=1):

		''' The data from the machine is processed in this function, inserted into a dictionary of dictionaries. It needs specification of which column all the data is. Default values have been given, based on our machine. 

		returns the dictionary of dictionaries (the data)
		and an array containing the well numbers, this is used in the iteration later'''


		my_data = genfromtxt(self.filename, delimiter=',')

		# !!! It currently reads in all the wells, but sometime in the future it would be nice to be able to tweak this, especially if you only loaded a fraction of the wells. 
		#wells=40 # user needs to be able to modify it


		self.array = []
		for i in xrange(1,wells):
			self.array.append('Well %s' % i)

		# creating a dictionary of dictionaries
		self.data = {}
		for x in xrange(1,wells): #USER INPUT = HOW MANY WELLS. 
			self.data['Well %s' % x] = {}		

		# reading in the relevant data and manipulating it 
		for i in xrange(1, my_data.shape[0]):	# ASSUMPTION: header is one line. 
			#print i, my_data[i,column_temp], round((my_data[i,column_temp]*1)+offset_temp, 2)
			my_data[i,column_temp] = (my_data[i,column_temp]*increment)+offset_temp
	
		
		# inserting into dictionary of dictionaries
		for values in xrange(1, my_data.shape[0]):	# iterating over the array. (my_data.shape[0] = nr of rows in the array)
			if my_data[values,columnwell] in xrange(1,wells): 
				 
				self.data['Well %s' % int(my_data[values,columnwell])][my_data[values,column_temp]]=my_data[values,yval]
				#sorted(self.data['Well %s' % int(my_data[values,columnwell])])

		self.emit(QtCore.SIGNAL('datadict'), self.data)
		return self.data, self.array


#################################################################################################################################
# ANALYSIS PART


	def slopetime(self):

		'''This function fits a sliding window of varying window size along the data, and evaluates where along the data we have a steepest increase by fitting 
		a line to each window and determining the corresponding r_value*slope value (so its normalized). 

		It starts with a window size 7, and then slides along the data with 1 increment and evaluates where we have the largest slope*r value. It then takes this 
		max value, and appends to a list, and simultaneously also makes a note of the window size. It does this for different windowsizes, with increment 2, up 
		until 35. We then get two lists, each windowsize has a corresponding r value which matches the max(slope*r) value. 

		ASSUMPTION: max(slope*r) is where we will find the melting temperature, as a result of the nature of the experiment. '''

		
		topr = []
		wsize = []
		for w_size in range(7,35,2): # We are trying different window sizes 
			# Clear these for each new window size 
			intercep = []
			slopes = []
			r = []
			wsize.append(w_size) # make note of which window size we are currently investigating
			for index in range(int((w_size-1)/2),int((len(yvalues)-(w_size-1)/2))): # And sliding these different window sizes along the array 
				ys = yvalues[index-((w_size-1)/2):index+((w_size-1)/2)] 
				xs = temperature[index-((w_size-1)/2):index+((w_size-1)/2)]

				slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xs,ys) # Calculating the slope at each window point

				# Appending to arrays 
				slopes.append(slope)
				intercep.append(intercept)
				r.append(r_value)
				multiply = [a*(b**2) for a,b in zip(slopes,r)]

			
			# by definition, where the melting temperature can be found, the slope*r will be greatest, so for all the windows of size w_size, make note
			# of the r_value which corresponds to the highest value of (slope*r). This will be where the melting temperature can be found. Also, before
			# we made  a note of the w_size. W_size and topr works as keys, to determine which windowsize will give the best melting temperature (using topr)
			topr.append(r[multiply.index(max(multiply))])

		return wsize, topr 

	def rtime(self):
		
		''' This function uses the topr and wsize arrays from slopetime. Empirically it has been determined that for a value of 
		r that roughly corresponds to 0.996 gives the best value for the melting temperature. However, this is an approximate number, so we need to find
		the windowsize that roughly correspondsself.canvas.draw() to an r value of 0.996, or at least the transition between 0.997 to something lower. The r value does not 
		drop linearly, i.e. it could be a drop down to 0.996 and then it would increase, and then later drop again. We want this "latest" drop, because
		this tends to be the most accurate. Hence, we extract all the values between 0.9954 and 0.9971 and then choose the windowsize that is the largest, 
		this will from empirical considerations be the best one for analysis.  '''


		window, r = self.slopetime()
		(a, b) = self.slopetime()
		c = zip(a,b)

		try:
			dough = [(a,b) for (a,b) in c if 0.9954<= b <= 0.9971] # Extract all rvalues with the corresponding wsize that have an rvalue between 0.9954 & 0.9971
			aa = [a for (a,b) in dough]
			#print "Value of r is:", max(aa)
			windowsize = max(aa) # From the previous extraction, choose the one with the biggest windowsize. 


		except ValueError: # means there is no value between 0.9954 and 0.9971. Just find the closest one.
			hello = min(enumerate(r), key=lambda x: abs(x[1]-0.996))
			#print "The value of r:", hello[1]
			windowsize = window[hello[0]]


		#print "optimal window size is then:", windowsize

		return windowsize 

	def windowtime(self):

		'''Optimal windowsize has been determined by the previous 2 functions, here we will implement the results from previous analysis, i.e. the 
		optimal windowsize and determine the melting temperature. The window function will "calculate the derivative" by taking subsequent and previous 
		indices and subtract them from the current y-value. To enhance the analysis, we then fill in the data using a quadratic spline (to avoid singularities)
		We define the melting temperature to be at the absolute peak of this spline curve. '''


		windowz = self.rtime() # obtain window size from dynamic window function
		lw = (windowz-1)/2
		rw = (windowz-1)/2

		# arrays that will contain the values of the window calculation 
		wwindow = []
		xxwindow = []

		# for a given window size, calculate the difference between each points 
		for index in range(int(lw),len(yvalues)-(int(rw)+1)):

			# do left and right window separartely. enables you to define an asymmetric window if needed
			leftfraction = []
			rightfraction = []

			for n in range(lw):
				leftfraction.append((yvalues[index]-yvalues[index-(n+1)]))

			for d in range(rw):
				rightfraction.append((yvalues[index]-yvalues[index+(d+1)]))

			wwindow.append((sum(leftfraction)-sum(rightfraction)))
			xxwindow.append(temperature[index])
			#print temperature[index],type(temperature[index])
		
		
		# SPLINE OF WINDOW 
		window_spl = scipy.interpolate.UnivariateSpline(xxwindow, wwindow, k=2)
		xwindow = linspace(min(xxwindow),max(xxwindow),1000)


		#max_array = scipy.signal.argrelmax(numpy.array(window_spl(xwindow)),order=15)
		#print max_array
		lists = list(window_spl(xwindow)) # make the spline array into a list
		maxis = lists.index(max(lists)) # find absolute maxima 
		#print "The melting temperature is suggested to be: ", xwindow[maxis]
		
		return xwindow, window_spl(xwindow), xwindow[maxis], maxis, max(window_spl(xwindow)) # returns spline wvalues, spline yvalues, the melting temperature, index of melting temperature, and the maximum of the spline yvalues for normalization later

	def peaktime(self):
		
		''' Sometimes, the denaturation is a 2-step process, this manifests as two peaks in the window-spline function. 
		To facilitate subsequent analysis, this function will detect the two peaks and return the temperature corresponding to the two peaks. '''

		x, vector, th, b, c = self.windowtime()
		scaling = len(x)/len(temperature)


		yheight = 0.02 # tweak this depening on window size?
		if self.rtime() > 11:
			windows = 6
		else:
			windows = 4

		maxima = []
		vector = list(vector/c)

		for yval in range(len(vector)):
			try:
				if vector[yval] > 0.3: # needs to be larger than 0.3, all other peaks are insignificant as small peaks 
					if vector[yval] > (vector[yval+(windows*scaling)]+yheight) and vector[yval] > (vector[yval-(windows*scaling)]+yheight):
						maxima.append(vector[yval])

				#elif 0 <= vector[yval] <= 0.3: # detecting unusual activity 

			except IndexError: # means we came to an end
				pass

		maxima.sort()
		good_indices = [ vector.index(fit) for fit in maxima ] # for the maxima, it finds the maximum index but also those around it. 
		groups = self.group_consecutives(sorted(good_indices), step=1)
		#print groups # for the sake of hewl-gi data

		#print "nr of peaks", len(groups)

		try:
			peak_1 = groups[0]

			pe_y1=[]
			pe_x1=[]
			for i in peak_1:
				pe_y1.append(vector[i])
				pe_x1.append(x[i])
			a= pe_x1[pe_y1.index(max(pe_y1))]

			if len(groups) > 1:
				peak_2 = groups[1]
				pe_y2=[]
				pe_x2=[]

				for i in peak_2:
					pe_y2.append(vector[i])
					pe_x2.append(x[i])

				f= pe_x2[pe_y2.index(max(pe_y2))]
			else:
				f = ''

		except ValueError: # means we have a peak really at the end 
			a = ''
			f = ''
		return a, f
		
	def group_consecutives(self, vals, step=1):
	    """Return list of consecutive lists of numbers from vals (number list)."""
	    run = []
	    result = [run]
	    expect = None
	    for v in vals:
	        if (v == expect) or (expect is None):
	            run.append(v)
	        else:
	            run = [v]
	            result.append(run)
	        expect = v + step
	    return result		

		
	def find_nearest(self,array,value):
		idx = (numpy.abs(array-value)).argmin()
		return array[idx]

	def peakanalysis(self, yarray, index_th, scaling, treshold):
		''' For analysing whether peaks identified are true peaks.'''

		try:
			a = index_th-(treshold*scaling)
			b = index_th+(treshold*scaling)

			if a<0: # means its so close to the edge it will go around
				print "Too close to the edge?"
				return "WARNING, Tm might not be true"

			else:
			# evaluation of the peak 
				before = yarray[index_th]-yarray[a]
				after = yarray[index_th]-yarray[b]
				
				if before > 0 and after >0: # true peak, values are smaller on either side
					return ''
				elif before <0 or after <0:
					print "Too close to the edge?"
					return "WARNING, Tm might not be true"

		except IndexError: # peak top is too close to the end of the array (high T), i think when b is a number longer than the array length
			print "NOOO"
			return "WARNING, Tm might not be true"

#################################################################################################################################

	def updategraph(self, status, *args):

		''' This is where it all  happens. The dictionary is iterated over (i.e. each well is processed) and the data-analysis functions are called, denaturation/two peak analysis is conducted and finally plotted. Also writes the resultsfile as it is iterating.

			as you can see, it calls the function that inserts the data into a dictionary. The parsefile function returns the dictonary of dictionaries (i.e. full data) and an array of all the wells, it is used the in iteration below. '''

		
		data, array = self.parsefile(column_temp=self.args[0]-1,offset_temp=self.args[1], yval=self.args[2]-1, columnwell=self.args[3]-1, wells=self.args[4]+1, increment=self.args[5]) # reading in the data, using the user comments into the correct datastructure (dictionary of dictionaries)

		#print data['Well %s' % 2]
		#############################################

		# RESULTS FILE, INITIATE AND WRITE FIRST ROW
		name = os.path.basename(self.filename)
		f = open("results.csv", "w")
		c = csv.writer(f)
		c.writerow(['Columntemperature', 'Offset_temperature', 'Filter','Columnwell'])
		c.writerow([self.args[0]-1,self.args[1],self.args[2]-1,self.args[3]-1])
		c.writerow(['Well Number', 'Calculated_Th', 'Peak_2/Warning','User_Input'])

		############################################
		#f = open(self.filename+'_results'+'.txt', 'a')
		#f.write('%s\t %s\t %s\t %s\n' % ('Well', 'Calculated_Th', 'Peak_2/Warning','User_Input'))


		for i in array:
			#################################################################################################
			# EXTRACTING DATA
			global yvalues, temperature # assignign global values 
			yvalues1 = data[i].values() #yvalues 
			temperature1 = data[i].keys() # Temperature 


			# IDK, BUT YOU NEED THIS IN ORDER TO GET THE DICTIONARY IN ORDER WHEN YOU HAVE INSERTED FLOATS AS A VALUE FOR INCREMENT OR OFFSET TEMPERATURE
			zipped = zip(temperature1, yvalues1)
			zipped.sort(key = lambda t: t[0])
			a = [list(t) for t in zip(*zipped)]
			temperature = a[0]
			yvalues = a[1]

			#print min(temperature), max(temperature)

			

			#################################################################################################
			# SPLINE
			spl = scipy.interpolate.UnivariateSpline(temperature, yvalues)
			x = linspace(min(temperature),max(temperature),1000)
	
			#################################################################################################
			# EXTRACING VALUES FROM ANALYSIS
			x, y, th, index_th , xth = self.windowtime()
			print i, ", The melting temperature is suggested to be: ", th
			#################################################################################################
			# NEED TO WARN USER WHEN THERE IS NO BASELINE (UNSTABLE PROTEIN)
			# either within a certain percentage of the end, or evaluate peak analysis, like done in peak_2!! But change treshold 
			# to 2 degrees

			self.comment = ''
	
			#################################################################################################
			# DENATURATION
			count = 0
			for points in y/xth:
				if points > 1.5 or points < (-1.0):
					count = count + 1

			#################################################################################################
			# TWO STAGE DECAY 
			
			peak_1, peak_2 = self.peaktime()

			# just make sure whatever peak that comes out of there (usually named Peak_1) is taken care of 
			
			if peak_1 != '':
				if (th-3) <= peak_1 <= (th+3):
					peak_1 = ''
		
			if peak_2 != '':
				if (th-3) <= peak_2 <= (th+3):
					peak_2 = ''

			########################################################################################################
			# FOR THE TABLE LATER

			if peak_2 != '' and peak_1 == '':
				#print tuple (sorted(round(th,1),round(peak_2,1)) )
				self.th_data.append( tuple (sorted( [round(th,1),round(peak_2,1)] ) ) )

			elif peak_1 != '' and peak_2 == '':
				#print tuple( sorted(round(th,1),round(peak_1,1)) )
				self.th_data.append( tuple( sorted( [round(th,1),round(peak_1,1)] ) ))

			elif count > 500:
				self.th_data.append('N/S')

			elif len(self.comment) > 1:
				self.th_data.append(((round(th,1),'ATTN')))

			else:
				self.th_data.append(round(th,1))
			#################################################################################################
			# PLOTTING AND WRITING DATA TO FILE

			####################################################

			# COLOURS
			scattercolor = [0,99.0/255.0,136.0/255.0]
			windowcolor = [159.0/255.0,161.0/255.0,97.0/255.0]
			derivativecolor = [126.0/255.0,49.0/255.0,123.0/255.0]
			peaks = [170.0/255.0,43.0/255.0,74.0/255.0]

			###################################################

			# MAIN PLOT
			ax = self.figure.add_subplot(1,1,1)
			ax.scatter(temperature, yvalues/(max(yvalues)), alpha=0.6, color=scattercolor)
			ax.plot(x, y/xth, alpha = 0.5, color=derivativecolor)
			axvline(th, color=peaks, alpha=0.7, linewidth=1.5) # MELTING TEMPERATURE()
			axvline(th - 0.5 * self.rtime(), color=windowcolor, alpha=0.7) # window begin
			axvline(th + 0.5 * self.rtime(), color=windowcolor, alpha=0.7) # window end
			txt = "T:"+str( '%1.1f' % th)+", "+self.comment

			# EXTRA PEAK DETECTION
			if peak_2 != '' and peak_1 == '':
				axvline(peak_2, color=peaks, alpha=0.7, linewidth=1.5)
				txt = "T:"+str(round(peak_2,1))+", "+str('%1.1f' % th)+self.comment
			elif peak_1 != '' and peak_2 == '':
				axvline(peak_1, color=peaks, alpha=0.7, linewidth=1.5)
				txt = "T: "+str(round(peak_1,1))+", "+str('%1.1f' % th)+self.comment

			# DENATURED
			if count > 500: # means more that 50% of data looks FANIIIII
				txt = 'NO SIGNAL FOUND'#, or maybe it is:'+str(round(th,2))		
			
			# AXES LABELS
			ylabel('Normalized Intensity')
			xlabel('Temperature ($^o$C)')
			text(min(temperature),1.1 , ("%s, " % i)+txt)
			ylim(-0.1,1.2)
			#ticklabel_format(style='sci', axis='y', scilimits=(0,0) 

			# show figure
			self.canvas.draw()

			#################################################################################################
			# SAVING THE GRAPH 
			if self.savfil.isChecked() == True:
				savefig('Well '+str(self.b[i])+'.png')
				self.canvas.figure.clf()
			else:
				self.canvas.figure.clf()
			#################################################################################################
			# AUTOMATIC MODE OR NOT	
				
			if status == 0:
				yield 
			#write to file only after they pressed continue! 
				c.writerow([i, round(th,1), str(peak_2)+self.comment, self.edit1.text()])
				self.edit1.setText('')
			#	continue

		#######################################################################################################################
		# LAST MINUTE FIX
		# wait, do I need this?
		self.canvas.figure.clf()
		f.close()

		# CALL TABLE FUNCTION TO DISPLAY TABLE, first we need to show the popup for the insertion of values
		self.widgetT.show() # display popup window


	def tempTable(self):

		''' This function first considers the suggestions by the user (the new temperature), creates a new dictionary based on this (again, the well is the key and the temperatures are the values). It then creates a colour coded table based on the reference temperature and what is considered a signficant change. 

		I think it has 5 different shades of change. It's a dynamic range, so it will adjust for each dataset (i.e. the darkest blue is always the largest increase) '''

		resultsfile = open("results.csv", 'r')
		datareader = csv.reader(resultsfile)
		next(datareader)

		count = 0
		# We will be modifying the self.meltingt data, instead of creating a whole new array. Userinputs will be replacing the program suggestions
		for row in datareader:
			#print row
			if count > 2:
					try:
						if len([ elem for elem in row[3].split(",") ]) > 1:
							print tuple([ float(elem) for elem in row[3].split(",") ])
							self.th_data[count] = tuple([ float(elem) for elem in row[3].split(",") ])
			
						else:
							self.th_data[count] = [ float(elem) for elem in row[3].split(",") ][0]
							print [ float(elem) for elem in row[3].split(",") ]

					except ValueError:
						pass

			count = count + 1

		self.emit(QtCore.SIGNAL('colordict'), self.th_data)
		# reference t and sig fig comes from the user input from the popup window


		#################################################################


		# OBTAINING REFERENCE AND SIG.CHANGE TEMPERATURE

		try: 
			ref = self.tableargs[0] # reference temperature, from userinput 
			minsigdelta = self.tableargs[1] # smallest significant change
		except: # this is the instance when you call the table from the menu instead, when re-reading in the results. Then it fetches the values for these parameters from an alterative window
			ref = float(self.referencet)
			minsigdelta = float(self.sigfigchange)


		################################################################

		#print "OK THIS IS THE TEMPTS THAT GO INTO THE TABLE", self.th_data


		# LABEL FOR COLUMNS AND ROWS
		labelc = ['1', '2', '3', '4', '5', '6', '7', '8', '9',
		          '10', '11', '12']
		labelr = ['',' A ', ' B ', ' C ', ' D ', ' E ', ' F ', ' G ', ' H ' ]

		cpool = [ (0.97,0.54,0.54),(0.83,0.25,0.25) , (0.70,0.08,0.08), (0.51,0.02,0.02), (0.43,0.05,0.05) ,
              (0.3,0.02,0.02), (0.0078,0.22,0.35), (0.015,0.352,0.55) ,(0.019,0.439,0.70) ,(0.211,0.56,0.75) , (0.65,0.66,0.811) ,
              (0.65,0.74,0.86) ,(0.65,0.74,0.86)]
		cmap3 = col.ListedColormap(cpool[0:12], 'indexed')
		cm.register_cmap(cmap=cmap3)



		# IF the user has only loaded <96 wells, we need to populate them so that we do not get an error in the table 
		poo = [1]*40
		print "the length is", len(self.th_data)
		print "the length of test is", len(poo)
		if len(self.th_data) < 96:
			96-len(self.th_data)
			for addintegers in range(1,97-len(self.th_data)):
				self.th_data.append('N/A')
		else:
			pass

		print len(self.th_data)









		# COLOUR TIME 
		colors = [[(0.95, 0.95, 0.95) for c in range(12)] for r in range(9)] 
		maxblue = max([x for x in self.th_data if type(x) in (numpy.float64, float)]) # ADJUST TEMPERATURE RANGE FOR MAXIMUM STABILISATION
		maxred = min([x for x in self.th_data if type(x) in (numpy.float64, float)]) # ADJUSTINGT TEMP RANGE FOR MAXIMUM DESTABILISATION
		change_b = [d for d in numpy.linspace(minsigdelta, int(maxblue-ref), 6)] # CREATE KEYS FOR THE DICTIONARY
		change_r = [d for d in numpy.linspace(minsigdelta, int(ref-maxred), 6)] # - || -
		red = [(0.99,0.99,0.70), (0.99,0.85,0.46),(0.99,0.698,0.29), (0.99,0.55,0.23), (0.94,0.23,0.13), (0.74,0,0.15)]
		blue = [(0.65,0.74,0.86),(0.65,0.66,0.811),(0.211,0.56,0.75),(0.019,0.439,0.70), (0.015,0.352,0.55), (0.0078,0.22,0.35)]
		cd_red = dict(zip(change_r, red)) # CREATE DICTIONARY
		cd_blue = dict(zip(change_b, blue)) # CREATE DICTIONARY
		cell_text = [['','','','','','','','','','','',''], self.th_data[0:12], self.th_data[12:24],self.th_data[24:36], self.th_data[36:48], self.th_data[48:60], self.th_data[60:72], self.th_data[72:84], self.th_data[84:96]]









		# COLOUR THE CELLS ACCORDING TO TEMPERATURE
		for lists in range(1,9,1): # skip empty first row
			for values in range(len(cell_text[lists])):
				try:
					th_change = float(cell_text[lists][values])-ref
					if th_change > minsigdelta: # the significant change should be able to be user defined 
						colors[lists][values] = cd_blue[self.find_nearest(numpy.asarray(change_b), th_change)] 
						#colors[lists][values] = (0.55, 0, 0)
					elif th_change < -minsigdelta: 
						colors[lists][values] = cd_red[self.find_nearest(numpy.asarray(change_r), numpy.abs(th_change))]
				except ValueError: # DENATURATION
					if len(list(cell_text[lists][values])) == 3 :
						colors[lists][values] = (0.5,0.5,0.5) # DENATURED
					else:
						cell_text[lists][values] = "\n".join(str(x) for x in cell_text[lists][values]) # SHoudl also work for warning message
						colors[lists][values] = (0.61,0.7,0.60)
				except TypeError: # TWO STAGE MELTING PROCESS
					#tuples = list(cell_text[lists][values]) # this is a tuple, so we need to un
					colors[lists][values] = (0.61,0.7,0.60)
					cell_text[lists][values] = "\n".join(str(x) for x in cell_text[lists][values])


		
		################################################


		# CREATE TABLE

		# title will be filename, uncomment if needed:
		tablename = os.path.basename(self.filename)
		#title(str(tablename)+" Screen")

		# TABLE PROPERTIES 
		lightgrn = (0.95, 0.95, 0.95)
		tab = table(cellText=cell_text,
		            rowLabels=labelr,
		            colLabels=labelc,
		            rowColours=[lightgrn]*16,
		            colColours=[lightgrn]*16,
		            cellColours=colors,
		            cellLoc='center',
		            loc='upper left')

		tab.scale(1,1.2)
		axis('off')
		table_props = tab.properties()
		table_cells = table_props['child_artists']
		for cell in table_cells: cell.set_height(0.1)

		#self.canvas.figure.clf()
		self.canvas.draw()

		#Saving the figure
		savefig(str(tablename)+"_Table"+'.png')


		# FORGOTTEN THEIR USED, TOO AFRAID TO DELETE
		#tab.auto_set_font_size(False)
		#tab.set_fontsize(16)
		#ax = self.figure.add_subplot(2,1,1)
		#colorbar.ColorbarBase(ax, cmap=cmap3)
		#ax.pcolormesh()
		#self.figure.set_size_inches(12,10)


		#######################################################################################################################
			
	def start(self):
		''' It starts the generator function.  '''
		self.widget.show() # the user should insert values when the start button is pressed 
		self._generator = self.updategraph(status=0)
		

	def continues(self, event):
		''' This is the generator function (I think), it forwards the iteration stepwise or automatically. '''
		# need to pass the value of the flick button to this thing 
		#print self.chick.isChecked()

		self.button2.setText("Continue")

		try: # here automatic mode is implemented, depending on the flick switch
			if self.chick.isChecked() == True:
				for value in self._generator:
					next(self._generator)
				#next(self._generator) # it starts from the very beginning 
			elif self.chick.isChecked() == False:
				next(self._generator)

		except StopIteration:
			print "it stopped"
			return

	def tableParameters(self, *args):
		self.tableargs = list(args) # now its a list
		return self.tableargs

	def about(self):
		''' If you want to insert details about the program, then add things here.'''

		QtGui.QMessageBox.about(self, "About",
"""Please cite our amazing paper Groftehauge, 2014, Acta D when using the program. """
)


class InputParametersPopup(QtGui.QWidget):
	''' This class contains the popup window where you insert all the files and gives the column values.'''
	def __init__(self):

		QtGui.QWidget.__init__(self)

		######################################


		self.setWindowTitle("Initial Settings")
		grid = QtGui.QGridLayout()
		self.setLayout(grid)
		self.setGeometry(1000, 300, 200, 400)
		self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
		self.m_widget = QtGui.QWidget(self)
		

		#####################################
		

		# FILE BUTTONS 

		self.rawdata = QtGui.QPushButton('Raw data')
		self.connect(self.rawdata, SIGNAL('clicked()'), self.sendingsignal1)
		grid.addWidget(self.rawdata, 1, 0)


		self.solutionsfile1 = QtGui.QPushButton('Solutes file (optional)')
		self.connect(self.solutionsfile1, SIGNAL('clicked()'), self.sendingsignal)
		grid.addWidget(self.solutionsfile1, 2,0)

		self.results1 = QtGui.QPushButton('Results from previous run (optional)')
		self.connect(self.results1, SIGNAL('clicked()'), self.sendingsignal2)
		grid.addWidget(self.results1,3,0)


		#######################################


		# CSV FILE COLUMNS

		# COLUMN TEMPERATURE 
		self.column = QtGui.QLineEdit()
		self.column.setPlaceholderText("1")
		#self.column.setText("0")
		grid.addWidget(self.column,4,1)
		lbl1 = QtGui.QLabel('Column containing temperature values ')
		lbl1.move(100,100)
		grid.addWidget(lbl1,4,0)

		# OFFSET TEMPERATURE
		self.offt = QtGui.QLineEdit()
		self.offt.setPlaceholderText("23")
		grid.addWidget(self.offt,5,1)
		lbl2 = QtGui.QLabel('Temperature offset ')
		lbl2.move(-3,10)
		grid.addWidget(lbl2,5,0)

		# INTENSITY DATA 
		self.int = QtGui.QLineEdit()
		self.int.setPlaceholderText("4")
		grid.addWidget(self.int,6,1)
		lbl3 = QtGui.QLabel('Column containing the intensity data ')
		# lbl3.move(16,10)
		grid.addWidget(lbl3,6,0)

		# WELL DATA
		self.well = QtGui.QLineEdit()
		self.well.setPlaceholderText("2")
		grid.addWidget(self.well,7,1)
		lbl4 = QtGui.QLabel('Column containing the well numbers')
		# lbl4.move(16,10)
		grid.addWidget(lbl4,7,0)

		# HOW MANY WELLS ARE LOADED
		self.wellnr = QtGui.QLineEdit()
		self.wellnr.setPlaceholderText("96")
		grid.addWidget(self.wellnr,8,1)
		lbl5 = QtGui.QLabel('Wells loaded')
		# lbl4.move(16,10)
		grid.addWidget(lbl5,8,0)

		# TEMPERATURE INCREMENTS
		self.tincr = QtGui.QLineEdit()
		self.tincr.setPlaceholderText("1")
		grid.addWidget(self.tincr,9,1)
		lbl6 = QtGui.QLabel('Temperature increments')
		# lbl4.move(16,10)
		grid.addWidget(lbl6,9,0)		


		#########################################

		# TWO BUTTONS

		# CONFIRM BUTTON
		self.buttonpress = QtGui.QPushButton('Apply')
		self.connect(self.buttonpress, SIGNAL('clicked()'), self.buttonClicked)
		grid.addWidget(self.buttonpress)

		# ClOSE BUTTON
		self.closebutton = QtGui.QPushButton('Close')
		self.connect(self.closebutton, SIGNAL('clicked()'), self.close)
		grid.addWidget(self.closebutton)


		#########################################


		# All these functions below are simply to recieve the filenames back from the main window to be displayed in the popupwindow. 

	def reciever_s(self, poo):
		#print "this is poo", poo
		self.solfilename = self.path_name(poo)

	def reciever_r(self, poo):
		#print "this is poo", poo
		self.resultsname = self.path_name(poo)

	def reciever_rawdata(self,poo):
		#print "this is poo3", poo
		self.rawdataname = self.path_name(poo)

	def sendingsignal(self):
		self.emit(QtCore.SIGNAL('opensolution'), "hello")
		self.solutionsfile1.setText(self.solfilename)

	def sendingsignal1(self):
		self.emit(QtCore.SIGNAL('openrawdata'), "hello")
		self.rawdata.setText(self.rawdataname)

	def sendingsignal2(self):
		self.emit(QtCore.SIGNAL('openresults'), "hello")
		self.results1.setText(self.resultsname)

	def path_name(self,path):
		''' Splits the pathname so that only the filename is shown in the window/on the buttons. '''

		head, tail = ntpath.split(path)
		return tail or ntpath.basename(head)

	def buttonClicked(self):
		''' Emits a signal carrying data: the column values etc. '''
		try:
			self.emit(QtCore.SIGNAL('mySig'), int(self.column.text()),float(self.offt.text()), int(self.int.text()), int(self.well.text()),int(self.wellnr.text()),float(self.tincr.text()) )
		except ValueError:
			QtGui.QMessageBox.critical(self, "Error", "Re-check your values.")


class TablePopup(QtGui.QWidget):
	''' This refers to the the little window that popups after the iterations asking for the signficant change in temperature and reference temperature. This is essentially the same as the one
	you can access in the menu (under data). Not sure what to do about this, I guess you can make the code more effective by calling this function. I have kept it for now, as it allows you 
	show the table and change the reference temperature around a bit. '''

	def __init__(self):
		QtGui.QWidget.__init__(self)
		###################################
		self.setWindowTitle("Creating a table")
		self.m_widget = QtGui.QWidget(self)
		dope = QtGui.QVBoxLayout(self.m_widget)
		####################################

		#  REFERENCE TEMPERATURE INPUT
		self.reft = QtGui.QLineEdit()
		self.reft.setPlaceholderText("0")
		dope.addWidget(self.reft)
		lbl1 = QtGui.QLabel('Reference temperature')
		lbl1.move(15,0)
		dope.addWidget(lbl1)

		####################################

		# SIGNFICANT TEMPERATURE 
		self.sigd = QtGui.QLineEdit()
		self.sigd.setPlaceholderText("1")
		#self.offt.setText("23")
		dope.addWidget(self.sigd)
		lbl2 = QtGui.QLabel('The smallest change in T that is considered significant')
		lbl2.move(16,10)
		dope.addWidget(lbl2)

		####################################


		# CONFIRM BUTTON
		self.buttonpress_table = QtGui.QPushButton('Apply')
		self.connect(self.buttonpress_table, SIGNAL('clicked()'), self.buttonClicked_table)
		dope.addWidget(self.buttonpress_table)

		# ClOSE BUTTON
		self.closebutton_table = QtGui.QPushButton('Close')
		self.connect(self.closebutton_table, SIGNAL('clicked()'), self.close)
		dope.addWidget(self.closebutton_table)

		# SHOW TABLE BUTTON
		self.showtablebutton = QtGui.QPushButton('Show table')
		self.connect(self.showtablebutton, SIGNAL('clicked()'), self.showTable)
		dope.addWidget(self.showtablebutton)

		#######################################

	def buttonClicked_table(self):
		''' Again, emits the signal with data to the main window class.'''
		try:
			self.emit(QtCore.SIGNAL('mySigT'), float(self.reft.text()),float(self.sigd.text()))
		except ValueError:
			QtGui.QMessageBox.critical(self, "C'mon, Retard", "Mate, insert some numbers")

	def showTable(self):
		self.emit(QtCore.SIGNAL('Table'))


class AnalysisPlotPopup(QtGui.QWidget):

	''' The class that deals with the plots (pH, temperature, waterfall plot). LIke the other classes it communicates with the mainwindow. '''

	def __init__(self):
		QtGui.QWidget.__init__(self)
		self.setWindowTitle("Plots for analysis")
		self.canvastime()
		self.marker = itertools.cycle(('x', '+', '^', 'o', '*', 'D', 's', 'p', 'H')) # this is for the plots later. 

	def canvastime(self, *data):
		self.m_widget = QtGui.QWidget(self)
		layout = QtGui.QVBoxLayout(self.m_widget)


		#####################################################################

		main = Dictionaries()
		try:
			#print list(data)[0]
			self.sdict = main.solutionsfile(list(data)[0])
		except:
			self.sdict = main.solutionsfile()



		a = self.sdict.values()
		solutions = []
		for i in range(len(a)):
			solutions.append(a[i][2])
		self.msg = set(solutions) # this is the solutes that have been clicked, i.e. the ones that the user wants to see
		self.remember = []

		self.plotdictionary_concentration = {}
		self.plotdictionary_pH = {}
		self.plotdictionary_cation = {}
		self.plotdictionary_conc_off = {}
		self.plotdictionary_ph_off = {}

		for i in self.msg:
			self.plotdictionary_concentration[i] = [], []
			self.plotdictionary_cation[i] = [],[]
			self.plotdictionary_pH[i] = [], []
			self.plotdictionary_conc_off[i] = [], []
			self.plotdictionary_ph_off[i] = [], []


		########## DROPDOWN BOX ##############################################
		self.button2 = QtGui.QToolButton()
		self.button2.setPopupMode(QtGui.QToolButton.MenuButtonPopup)
		self.button2.setText('Solutes')
		self.button2.setMenu(QtGui.QMenu(self.button2))
		layout.addWidget(self.button2)

		self.linear = QtGui.QRadioButton('linear')
		self.linear.setChecked(True)
		layout.addWidget(self.linear)

		self.log = QtGui.QRadioButton('log')
		layout.addWidget(self.log)

		##### SOLUTION BUTTONS ################################################

		for i in self.msg:
			self.button_i = QtGui.QPushButton(i)
			action_i = QtGui.QWidgetAction(self.button2)
			self.button_i.clicked.connect(partial(self.sorting, i) ) # 3D GRAPH
			self.button_i.clicked.connect(partial(self.cleargraph, i) ) # TH VS CONCENTRATION 
			self.button_i.clicked.connect(partial(self.correctedph, i) ) # TH VS. CORRECTED PH
			action_i.setDefaultWidget(self.button_i)
			self.button2.menu().addAction(action_i)


		##########  CANVAS ####################################################
		self.figure = plt.figure(figsize= (10,6), dpi=100)
		self.canvas =  mpl.backends.backend_qt4agg.FigureCanvasQTAgg(self.figure)
		layout.addWidget(self.canvas)
		self.canvas.setParent(self.m_widget) # main_widget parent instead of main_frame
		self.canvas.setFocusPolicy(Qt.ClickFocus)
		self.canvas.setFocus()
		self.canvas.draw()

		self.mpl_toolbar = NavigationToolbar(self.canvas, self)
		layout.addWidget(self.mpl_toolbar)

		#######################################################################


		# BUTTONS

		self.button_3D = QtGui.QPushButton('Waterfall Plot')
		self.button_3D.clicked.connect(self.graph3D)
		layout.addWidget(self.button_3D)

		self.button_conc = QtGui.QPushButton('Plot T vs. Concentration')
		self.button_conc.clicked.connect(self.graphconcentration)
		layout.addWidget(self.button_conc)

		#self.button_clear = QtGui.QPushButton('Cations')
		#self.button_clear.clicked.connect(self.plotcationtime)
		#layout.addWidget(self.button_clear)

		self.button_ph =  QtGui.QPushButton('Plot T vs. Corrected pH')
		self.button_ph.clicked.connect(self.graphpH)
		layout.addWidget(self.button_ph)

		self.button_clear = QtGui.QPushButton('Clear')
		self.button_clear.clicked.connect(self.clear)
		layout.addWidget(self.button_clear)

		#########################################################################

	def sorting(self, value):
		
		''' Creates a new dictionary of those wells that contain the data for a specific solution (when you want to plot all the wells containing Nickel) '''

		#self.remember.append(value)

		### EXTRACTION OF THE RELEVANT WELLS ###############################################################
		self.concentrationforplot = []
		self.wellsforplot = []

		a = self.sdict.values()
		b = self.sdict.keys()
		for i in range(len(a)):
			if a[i][2] == value:
				self.concentrationforplot.append(a[i][3])
				self.wellsforplot.append(b[i])

		self.verts = []
		for i in self.wellsforplot:
			self.data[i] # these dictionaries have the x values as the key and the intensities as the yvalues
			x_i = self.data[i].keys()
			z_i = self.data[i].values()/max(self.data[i].values())
			z_i[0],z_i[-1] = 0, 0
			self.verts.append(zip(x_i, z_i))
			y_i = [self.concentrationforplot[self.wellsforplot.index(i)]]

	def graph3D(self):

		''' Colours and plots the 3D waterfall plot.'''

		######## COLOURING OF THE RELEVANT WELLS ###########################################################
		self.figure.clf()
		ax = self.figure.gca(projection='3d')
		
		#print self.referencet
		#print self.meltingt
		#print self.diff

		try: 
			###############################################################################################################

			# COLOUR TIME 
			maxblue = max([x for x in self.meltingt if type(x) in (numpy.float64, float)]) # ADJUST TEMPERATURE RANGE FOR MAXIMUM STABILISATION
			maxred = min([x for x in self.meltingt if type(x) in (numpy.float64, float)]) # ADJUSTINGT TEMP RANGE FOR MAXIMUM DESTABILISATION
			change_b = [d for d in numpy.linspace(self.diff, int(maxblue-self.referencet), 6)] # CREATE KEYS FOR THE DICTIONARY
			change_r = [d for d in numpy.linspace(self.diff, int(self.referencet-maxred), 6)] # - || -
			red = [(0.97,0.54,0.54, 0.4), (0.83,0.25,0.25, 0.4),(0.70,0.08,0.08, 0.4), (0.51,0.02,0.02, 0.4), (0.43,0.05,0.05, 0.4), (0.3,0.02,0.02, 0.4)]
			blue = [(0.65,0.74,0.86, 0.4),(0.65,0.66,0.811, 0.4),(0.211,0.56,0.75, 0.4),(0.019,0.439,0.70, 0.4), (0.015,0.352,0.55, 0.4), (0.0078,0.22,0.35, 0.4)]
			cd_red = dict(zip(change_r, red)) # CREATE DICTIONARY
			cd_blue = dict(zip(change_b, blue)) # CREATE DICTIONARY

			###############################################################################################################

			zre = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] # this is the yvalues going straight up, for plotting the vertical lines (melting temperatures) 
			lines = []
			colors = []


			for i in self.wellsforplot:
				idx = int(i[5:])-1

				if self.meltingt[idx] == 'N/S': 
					colors.append((0.5,0.5,0.5, 0.4))
					th_x = [0] * 10

				elif type(self.meltingt[idx]) == tuple:
					if type(self.meltingt[idx][1]) not in (numpy.float64, float): # ATTN
						th_x = [self.meltingt[idx][0]] * 10


						th_change = float(self.meltingt[idx][0])-self.referencet

						if th_change > self.diff:
							colors.append(cd_blue[self.find_nearest(numpy.asarray(change_b), th_change)])

						elif th_change < -self.diff:
							colors.append(cd_red[self.find_nearest(numpy.asarray(change_r), numpy.abs(th_change))])

						else:
							colors.append((0.56,0.45,0.64, 0.4))

					

					else: # two stage
						colors.append((0.61,0.7,0.60, 0.4))
						th_x = [self.meltingt[idx][0]] * 10 # taking the smaller temperature, to highlight the dataset as funny/wrong 

				else:
					th_change = float(self.meltingt[idx])-self.referencet
					if th_change > self.diff:
						colors.append(cd_blue[self.find_nearest(numpy.asarray(change_b), th_change)])
					elif th_change < -self.diff:
						colors.append(cd_red[self.find_nearest(numpy.asarray(change_r), numpy.abs(th_change))])
					else:
						colors.append((0.56,0.45,0.64, 0.4))

					th_x = [self.meltingt[idx]] * 10


				lines.append(list(zip(th_x, zre)))

			poly = PolyCollection(self.verts, facecolor=colors, closed=False, edgecolor=None)
			line = LineCollection(lines, alpha=0.8, color='r')

		except AttributeError: # when there has been no iteration, user wants just to see the shape of the data, no colour no th
			poly = PolyCollection(self.verts,closed=False, edgecolor=None)

		#########################################################################################################################################################
		# PLOTTING ON A LINEAR OR LOG SCALE, DEPENDING ON THE STATUS OF THE RADIOBUTTONS


		if self.linear.isChecked()==True:
			ax.add_collection3d(poly, zs=self.concentrationforplot, zdir='y')
			ax.set_ylabel('Concentration (mM)')
			ax.set_ylim3d(max(self.concentrationforplot), min(self.concentrationforplot))
			try:
				ax.add_collection3d(line, zs=c, zdir='y') # adding the melting temperatures 
			except:
				pass
		elif self.log.isChecked()==True:
			c = [math.log(i) for i in self.concentrationforplot]
			ax.add_collection3d(poly, zs=c, zdir='y')
			ax.set_ylabel('Concentration (Log$_{10}$mM)')
			ax.set_ylim3d(max(c), min(c))
			try:
				ax.add_collection3d(line, zs=c, zdir='y')
			except:
				pass

		else:
			ax.add_collection3d(poly, zs=self.concentrationforplot, zdir='y')
			ax.set_ylabel('Concentration (mM)')
			ax.set_ylim3d(max(self.concentrationforplot), min(self.concentrationforplot))
			try:
				ax.add_collection3d(line, zs=c, zdir='y')
			except:
				pass
		#####################################################################################################


		# STUFF FOR THE FIGURE 
		#ax.add_collection3d(poly, zs=c, zdir='y')
		poly.set_alpha(0.4)
		ax.set_xlabel('Temperature ($^o$C)')
		ax.set_xlim3d(23, 100)
		ax.set_zlabel('Normalized  Intensity')
		ax.set_zlim3d(0, 1)
		self.canvas.draw()

		###################################################################################################

	###################################################################################################

		# RECIEVING DATA FROM THE OTHER CLASS, HOW SHOULD I MAKE THIS BETTER??? 
	def graphdata(self, data):
		self.data = data 
		
	def colouring(self, data): # this is the melting temperatures, we need
		self.meltingt = data
		print self.meltingt

	def values(self, ref, diff):
		self.referencet = ref
		self.diff = diff

	def cationtime(self):
		''' Function that will plot the ionic strength. Still at the rudimentary level.'''

		if len(set(self.remember)) == 0:
			QtGui.QMessageBox.critical(self, "Error", "You haven't clicked on any solutes")


		for value in set(self.remember):
			# value is the solute 
			print value
			self.plotdictionary_cation[value] = [], []
			a = self.sdict.values() # read in solution file

			for i in range(1,96):
				
				if self.sdict['Well %s' % i][2] == value:
					
					if a[i][6] != 0 and a[i][7] != 0: # check the cationic strength
						concentration = float(a[i][3])*int(abs(a[i][7])) # given that the cation has smaller charge 
						self.plotdictionary_cation[value][1].append(concentration) # HERE IS THE CATIONIC STRENGTH

						# need to go to keys and find out which well it corresp
						if type(self.meltingt[i-1]) not in (numpy.float64, float):
							if len(list(self.meltingt[i-1])) == 3:
								self.plotdictionary_cation[value][0].append(0)

						if type(self.meltingt[i-1]) == tuple:
							if type(self.meltingt[i-1][1]) not in (numpy.float64, float): # means ATTN
								self.plotdictionary_cation[value][0].append(self.meltingt[i-1][0])
							else: # means two stage
								self.plotdictionary_cation[value][1].append(self.meltingt[i-1][0]) # appending only the first tvalue 
						else:
							self.plotdictionary_cation[value][0].append(self.meltingt[i-1])

			print self.plotdictionary_cation[value][0], self.plotdictionary_cation[value][1]


	def plotcationtime(self):
		self.canvas.figure.clf()
		ax = self.figure.add_subplot(111)
		ax.set_ylabel('T$_h$ ($^o$C)')
		self.cationtime()
		for solute in set(self.remember):
			if solute != False:
				if len(self.plotdictionary_cation[solute][1]) > 4:
					ax.plot(ax,self.plotdictionary_cation[solute][1],self.plotdictionary_cation[solute][0], label=str(solute), marker='o')
				else:
					ax.scatter(ax,self.plotdictionary_cation[solute][1],self.plotdictionary_cation[solute][0], label=str(solute), marker='s', s=60)


		legend_font_props = FontProperties()
		legend_font_props.set_size('small')
		ax.set_xlabel('Cation Concentration (mM)')
		ax.legend(ax,loc='best', prop=legend_font_props)
		self.figure.canvas.draw()
		
						

	def cleargraph(self,value):
		self.remember.append(value)
		print self.remember
		return self.remember

	def clear(self):
		self.remember = []
	###################################################################################################

	def logdetection(self, array):

		''' Not sure if useful, but I have kept the code anyway. Detects logs given the input of an array. It's really buggy though '''
		s = sorted(array, reverse=True)
		d = [ round(s[i]/s[i+1]) for i in range(len(s)-1)]
		try:
			average = int(numpy.mean(numpy.asarray(d)))
		except ValueError: # this happens when no concentration has been given
			QtGui.QMessageBox.critical(self, "Error", "No supplied concentrations, unable to create 3D plot")

		g = Counter(d).most_common(1)
		logvalue = int(g[0][0])
		occurences = int(g[0][1])

		if occurences < 0.4 * len(d) or logvalue==1:
			return logvalue == 0
		else:
			return logvalue

	def find_nearest(self,array,value):
		idx = (numpy.abs(array-value)).argmin()
		return array[idx]

	def concentrationgraphs(self):

		''' Tries to extract the melting temperature and plot it against concentration. First pulls out all those temperatures that corresponds with the given solute you are wanting to plot. 
		Because the melting temperature can also be "DENATURED", "ATTN" or any other user comment the code is quite convoluted. For denatured solutes I just null the temperature. '''

		for value in self.msg:

			# create for loop, create datasets for the one we want to look at. 
			self.plotdictionary_concentration[value] = [], [] # think I am nulling this before. otherwise when you call it multiple times, the graph looks weird
			self.plotdictionary_conc_off[value] = [], []
			for i in range(1,97):
				if self.sdict['Well %s' % i][2] == value: # pulling out those values the user wants to plot 

					print self.meltingt[i-1], type(self.meltingt[i-1])

					# DENATURATION ##############
					if type(self.meltingt[i-1]) not in (numpy.float64, float): 
						if len(list(self.meltingt[i-1])) == 3:
							self.plotdictionary_concentration[value][1].append(0) # JUST APPENDING 0 T for now, what should I do?????
							self.plotdictionary_concentration[value][0].append(self.sdict['Well %s' % i][3])

							self.plotdictionary_conc_off[value][1].append(0)
							self.plotdictionary_conc_off[value][0].append(self.sdict['Well %s' % i][3])
					##############################
					
					# ATTN/TWO STAGE
					if type(self.meltingt[i-1]) == tuple: # means two stage or ATTN # THIS COULD POTENTIALLY BE WHERE IT GOES WRONG 

						if type(self.meltingt[i-1][1]) not in (numpy.float64, float): # means ATTN
							self.plotdictionary_concentration[value][1].append(self.meltingt[i-1][0])
							self.plotdictionary_concentration[value][0].append(self.sdict['Well %s' % i][3])
							print "ATTN"

							self.plotdictionary_conc_off[value][1].append(self.meltingt[i-1][0]) # append the melting temperature 
							self.plotdictionary_conc_off[value][0].append(self.sdict['Well %s' % i][3]) # the concentration


						else: # means two stage 
							self.plotdictionary_concentration[value][1].append(self.meltingt[i-1][0]) # appending only the first tvalue 
							self.plotdictionary_concentration[value][0].append(self.sdict['Well %s' % i][3])
							print "two stage"

							self.plotdictionary_conc_off[value][1].append(self.meltingt[i-1][1]) # append the other two stage denaturation value 
							self.plotdictionary_conc_off[value][0].append(self.sdict['Well %s' % i][3]) # the concentration

					###############################
					# NORMAL
					else:
						self.plotdictionary_concentration[value][1].append(self.meltingt[i-1])
						self.plotdictionary_concentration[value][0].append(self.sdict['Well %s' % i][3])

		
	def correctedph(self): # READS IN ALL VALUES, HENCE WILL BE PROBLEM IF ANYTHING IS NOT CORRECTLY INSERTED. 
		''' This is the function that extracts the pH value and the temperature that corresponds to the solute that the user wants to plot.
		'''


		for value in set(self.remember):
			if value != False:

				self.plotdictionary_pH[value] = [], []
				for i in range(1,97):
					if self.sdict['Well %s' % i][2] == value:


						# THREE SCENARIOS: DENATURED, ATTN/TWOSTAGE, NORMAL 

						# DENATURED ##########
						if type(self.meltingt[i-1]) not in (numpy.float64, float):
							if len(list(self.meltingt[i-1])) == 3:
								self.plotdictionary_concentration[value][1].append(0)
						######################

						else:
							######### TWO STAGE/ATTN ##############
							if type(self.meltingt[i-1]) == tuple:
								if type(self.meltingt[i-1][1]) not in (numpy.float64, float): # means ATTN
									th = self.meltingt[i-1][0]
								else: # two stage denaturation 
									th = self.meltingt[i-1][0]# ,-- not sure if this one works 
							######################################


							### NORMAL #############################
							else: 
								th = self.meltingt[i-1]

							print self.sdict['Well %s' % i][1], self.sdict['Well %s' % i][0]
							corrected = ((th-20)*float(self.sdict['Well %s' % i][1]))+float(self.sdict['Well %s' % i][0])
							self.plotdictionary_pH[value][0].append(corrected)
							self.plotdictionary_pH[value][1].append(th)
							########################################


	def graphconcentration(self):
		''' Function that plots melting temperature vs. concentration. If anything goes wrong with these plots it is usually here.'''
		

		print self.referencet

		self.canvas.figure.clf()
		ax = self.figure.add_subplot(111)
		#box = ax.get_position()
		#ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9 ])
		ax.set_ylabel('T$_h$ ($^o$C)')
		self.concentrationgraphs() # when you call this multiple times, is affected?
		
		for solute in set(self.remember):
			if solute != False: 

				# The following three lines are just to sort the plotting so they are plotted in order and it doesnt look funny when you try to plot it
				print (self.plotdictionary_concentration[solute][0],self.plotdictionary_concentration[solute][1])
				zipped = zip(self.plotdictionary_concentration[solute][0],self.plotdictionary_concentration[solute][1])
				zipped.sort(key = lambda t: t[0])
				a = [list(t) for t in zip(*zipped)]

				ax.plot(a[0],a[1], label=str(solute), marker=self.marker.next())

		if self.linear.isChecked()==True:
			ax.set_xlabel('Concentration (mM)')
		elif self.log.isChecked()==True:
			ax.set_xscale('log')
			ax.set_xlabel('Concentration (Log$_{10}$mM)')
		else:
			pass

		axhline(self.referencet, alpha=0.4)
		legend_font_props = FontProperties()
		legend_font_props.set_size('small')
		ax.legend(loc='best', prop=legend_font_props)
		self.figure.canvas.draw()



	def graphpH(self):
		''' Function that plots melting tempetature vs concencentration. Less buggy than the above one. '''

		self.canvas.figure.clf()
		ax = self.figure.add_subplot(111)
		ax.set_ylabel('T$_h$ ($^o$C)')
		self.correctedph()

		for solute in set(self.remember):
			if solute != False:
				# The following three lines are just to sort the plotting so they are plotted in order and it doesnt look funny when you try to plot it
				zipped = zip(self.plotdictionary_pH[solute][0],self.plotdictionary_pH[solute][1])
				zipped.sort(key = lambda t: t[0])
				a = [list(t) for t in zip(*zipped)]

				ax.plot(a[0],a[1], label=str(solute), marker=self.marker.next())
				#ax.plot(self.plotdictionary_pH[solute][0],self.plotdictionary_pH[solute][1], label=str(solute), marker=self.marker.next())
				#previous plotting, before we sorted it
		legend_font_props = FontProperties()
		legend_font_props.set_size('small')
		ax.set_xlabel('Corrected pH')
		ax.legend(loc='best', prop=legend_font_props)
		self.figure.canvas.draw()
		

path_to_namipy = os.path.dirname(os.path.abspath(__file__))
qApp = QtGui.QApplication(sys.argv)
app_icon = QtGui.QIcon()
app_icon.addFile(str(path_to_namipy)+'/fox.png', QtCore.QSize(256,256))
qApp.setWindowIcon(QtGui.QIcon(app_icon))
aw = ApplicationWindow()
aw.setWindowTitle("%s" % "NAMI")
aw.show()
sys.exit(qApp.exec_())
qApp.exec_()

