#This notebook will take a fits file and plot it using wcs coordinates. should look neat
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io.votable import parse
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 18})

def gauss(mx,mu,sigma,A):
    return A*np.exp(-(mx-mu)**2/2/sigma**2)

def fitsplot(image):
#file = '/home/gonzalo/astronomy_practice/Lynx_E_3.fits' #path of file where HST data was downloaded
	fitsfile = fits.open(image) #open the fits file
	wcs = WCS(fitsfile[0].header,naxis=2) #open wcs solutions
	data = fitsfile[0].data #get the data at each pixel

	#Plot a histogram
	y,x,_= plt.hist(data.flatten(),10000,alpha=1,label='data')
	#mean = np.mean(data.flatten())
	mx = (x[1:]+x[:-1])/2 # for len(x)==len(y)
	expected = (100,10,100) #guess parameters for gaussian : mean, std, amp

	#solve for parameters using scipy curvy fit, guess values and the data
	params,cov = curve_fit(gauss,mx,y,expected)
	sigma = np.sqrt(np.diag(cov))
	plt.plot(mx,gauss(mx,*params),color='red',lw=3,label='model') #plotting the double gaussian
	plt.legend()
	plt.xlim(0,20000)
	#print(params,'\n',sigma)

	#lower and upper limit for colorbar
	vl = params[0] - ( 3 * params[1]) 
	vm = params[0] + (3 * params[1])
	
	print(vl,vm)

	fig, axes = plt.subplots(figsize=(16,14),nrows=1, ncols=2, sharex=False,sharey=True,subplot_kw={'projection': wcs})

	#plotting parameters
	axes[0].imshow(data,cmap='gray') #plot it
	axes[0].set_title('M33')
	axes[0].set_xlabel('Right Ascention')
	axes[0].set_ylabel('Declination')

	axes[1].imshow(data[:,:],cmap='gray',vmin=vl,vmax=vm) #plot it
	axes[1].set_title('M33 w/ color adjustment')
	axes[1].set_xlabel('Right Ascention')

	plt.savefig('/home/gonzalo/astronomy_practice/test.png')
	plt.show()

print('What is the fits file path?')
image = input()
fitsplot(image)
