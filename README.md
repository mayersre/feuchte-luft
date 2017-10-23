feuchte-luft
============

There is an english version now and java code too, see my repository HumidAir

Feuchte Luft in Excel und python. Humid Air in Excel and python (most description in german) 

Excel version provides functions in excel for calculation of humid air readable in VBA

A port to Libre office should be easy as all is plain VBA

Python version is a library with the same functions for humid air calculations

All functions are from common Thermodynamics books or from lessons on humid air.

Humid air is treated like an ideal gas.

hair = hUMID air
(c) 2010 2011 2012 Reiner Mayers
This code comes under the GNU Copyleft (GPL) V3.0

Units with one exception :
  - all pressures in bar absolute
  - all Temperatures in degrees Celsius
The little letters behind the functions tell about the inputs :

   p       Pressure in bar absolute
   t       Temperature in degrees Celsius
   x       HumidityRatio in kg water per kg dry air
   h       Enthalpy in kJ/kg
   phi     Relative Humidity as a number 0.5 = 50% ; 1 = 100%
   pw      partial pressure of water
   pws     pw but saturated
   tDew    Dew point temperature in Celsius
   twb     Wet Bulb Temperature

