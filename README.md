# Playing with XDS predictions (using Numpy)
This is a tutorial at [ECACOMSIG computing school (Freudenstadt, August 2016)](http://www.mrc-lmb.cam.ac.uk/harry/ecacomsig/freudenstadt.html)

In this tutorial we will implement calculation of spot positions in 3D space (x,y,&phi;) from [XPARM.XDS](http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html#XPARM.XDS).

Reference paper: [W. Kabsch "Integration, scaling, space-group assignment and post-refinement" Acta Cryst. (2010). D66, 133-144](http://dx.doi.org/10.1107/S0907444909047374)

Test data: [Thaumatin / Diamond Light Source I04 user training](https://zenodo.org/record/10271#.V6vChZOLRBz)

If you want to run the script, you need to install
* Python2.7 with [CCTBX](http://cctbx.sourceforge.net/) and [Numpy](http://www.numpy.org/) (1.8 or higher)
   * You can use phenix.python if you installed [phenix-1.10.1](https://www.phenix-online.org/)
* [XDS and generate_XDS.INP](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Installation)
* [Adxv](http://www.scripps.edu/tainer/arvai/adxv.html)

## How to run the script
1. Get sample data
<pre>
mkdir th_8_2
cd th_8_2
wget https://zenodo.org/record/10271/files/th_8_2.tar.bz2
tar xvf th_8_2.tar.bz2
</pre>
2. Prepare XDS.INP and run XDS (Alternatively, you can use XPARM.XDS in this repository)
<pre>
mkdir xds
cd xds
generate_XDS.INP "../th_8_2_0???.cbf"
vi XDS.INP # Edit SPOT_RANGE= 1 27   JOB= XYCORR INIT COLSPOT IDXREF
xds_par
</pre>
3. Run the script
<pre>
phenix.python prediction.py ./XPARM.XDS 1.5 0.06928 1
# Usage:      prediction.py [XPARM.XDS] [d_min] [mosaicity] [frame numbers..]
</pre>
4. Use Adxv to see predictions with image
<pre>
adxv ../th_8_2_0001.cbf prediction_000001.adx
</pre>
