*This project is a fork of thermistor.sourceforge.net that I created in order to add features, data and to make the installation easier.*
# 1 Abstract

The project offers support for NTC thermistor calculations. The Steinhart-Hart equation is a mathematical model for these thermistors that seems to fit for a wide range of temperatures with high precision. Software to calculate the characteristic Steinhart-Hart coefficients based on temperature-resistance tables for given thermistors as well as functions allowing conversion of temperature values to resistance and vice versa is provided.

# 2 Description

A model for the resistivity of a semiconductor as a function of the temperature was found by Steinhart and Hart 1968 ([1]). The Steinhart-Hart law describes the absolute temperature T (in Kelvins) as a function of the NTC thermistor's resistivity (in Ω) according to the formula

<table align="center">
  <caption align="bottom">
  <em>    Steinhart-Hart polynom
  </em>
  </caption>
  <tbody><tr>
    <td><div align="center"><em><strong><sup>1</sup>/<sub>T</sub> = a<sub>0</sub> + a<sub>1</sub> · ln r + a<sub>3</sub> · (ln r)<sup>3</sup> </strong></em></div></td>
  </tr>
</tbody></table>

The figure below shows the typical graph of an NTC thermistors characteristic, giving the reciprocal temperature (in 1 / K) over the natural logarithm of the resistance (in Ω).

-> ![Characteristic curve of an NTC thermistor](https://raw.githubusercontent.com/epsilonrt/ntc/master/doc/html/ntcthermistor.png) <-

The constants a0, a1 and a3, also called Steinhart-Hart coefficients, vary depending on the type of thermistor. To support developer when creating temperature measurement applications, thermistor manufacturer often supply these constants for their products. They also publicate tables where resistivity of thermistor products for a wider range of temperature values are listed.

This project provides software to

* calculate temperature value for a given resistance of an NTC thermistor with given Steinhart-Hart coefficients,  
* calculate resistance value for a given temperature for an NTC thermistor with given Steinhart-Hart coefficients and  
* evaluate Steinhart-Hart coefficients for an NTC thermistor descibed by a temperature-resistance table.

Apart from the standard Steinhart-Hart equation other forms have been found. For application with lower CPU power a simplified form of the Steinhart-Hart equation can be used. 

<table width="241" border="0" align="center">
	<caption align="bottom">
	<em>  Simplified Steinhart-Hart polynom
	</em>
	</caption>
	<tbody><tr>
		<td width="235"><div align="center"><strong><em><sup>1</sup>&#8260;<sub>T</sub> = a<sub>0</sub> + a<sub>1</sub> · ln r</em></strong></div></td>
	</tr>
</tbody></table>

On the other hand a quadratic term can be inserted into the formula to increase accuracy giving the extended Steinhart-Hart equation

<table border="0" align="center">
	<caption align="bottom">
	<em> Extended Steinhart-Hart polynom </em>
	</caption>
	<tbody><tr>
		<td width="294"><div align="center"><strong><span class="Stil1"><sup>1</sup>/<sub>T</sub> = a<sub>0</sub> + a<sub>1</sub> · ln r + a<sub>2</sub>· (ln r)<sup>2</sup> + a<sub>3</sub> · (ln r)<sup>3</sup></span></strong></div></td>
	</tr>
</tbody></table>


An introduction to thermistors and the Steinhart-Hart polynom can be found at Wikipedia [[2](http://en.wikipedia.org/wiki/Thermistor)]. 

# 3 Software

Software is provided by this project for the calculations given in section algorithms. 

The classes / modules can be used for:

* calculation of Steinhart-Hart coefficients for an NTC thermistor whose characteristics is given as T-R table (utils/coeff).
* conversion from resistance to temperature (test/r2t),
* conversion from temperature to resistance (test/t2r),

The latter can be done for standard, simplified or extended Steinhart-Hart polynom.

## Installation

    git clone https://github.com/epsilonrt/ntc.git
    cd ntc
    make
    sudo make install

## Calculation of Steinhart-Hart coefficients

    cd ntc-data
    ntc-coeff murata-nxft15-10k.csv 
    Thermistor library version 1.0
    Copyright (C) 2007, 2013 - SoftQuadrat GmbH, Germany

    Steinhart-Hart coefficients
    a[0] = 9.310296797541951e-04
    a[1] = 2.308343095769287e-04
    a[2] = 3.001370069362199e-06
    a[3] = 5.407975166655454e-08


# 4 License

The software of this project is published under the GNU Lesser General Public License.

# 5 References

[1] J.S. Steinhart and S.R. Hart, “Calibration Curves for Thermistors,” Deep Sea Research, Vol. 15, 497-503, 1968.

[2] [Thermistors, NTCs and Steinhart-Hart equation on Wikipedia](http://en.wikipedia.org/wiki/Thermistor)
