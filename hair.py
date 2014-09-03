#
# hair = hUMID air
# (c) 2010 2011 2012 Reiner Mayers
# This code comes under the GNU Copyleft (GPL) V3.0
#
# Units with one exception :
#   - all pressures in bar absolute
#   - all Temperatures in degrees Celsius
#
# The little letters behind the functions tell about the inputs :
#
#   p       Pressure in bar absolute
#   t       Temperature in degrees Celsius
#   x       HumidityRatio in kg water per kg dry air
#   h       Enthalpy in kJ/kg
#   phi     Relative Humidity as a number 0.5 = 50% ; 1 = 100%
#   pw      partial pressure of water
#   pws     pw but saturated
#   tDew    Dew point temperature in Celsius
#   twb     Wet Bulb Temperature
#
# Constant Values
#
import math
#
import numpy as np
#
from scipy.optimize import brentq,newton
#
# R of Air
def RAir():
    """
    Gas constant of Air
    """
    _ret = None
    _ret = 287.05
    return _ret

#R of Steam
def RSteam():
    """
    Gas constant of Steam 
    """
    _ret = None
    _ret = 461.52
    return _ret

# RAir/RSteam = 0.622
def RA_RS():
    _ret = None
    _ret = RAir() / RSteam()
    return _ret

# Isobaric heat capacity of Air
def cpA():
    _ret = None
    _ret = 1.0046
    return _ret

# Evaporation enthalpy of water at 273.16K
def R0():
    _ret = None
    _ret = 2500.5
    return _ret

# Isobaric heat capacity of Steam
def cpS():
    _ret = None
    _ret = 1.88
    return _ret

# Heat capacity of water
def cw():
    _ret = None
    _ret = 4.192
    return _ret

# Enthalpy required to freeze water to ice 
def cE():
    _ret = None
    _ret = 2.07
    return _ret

# Heat capacity of ice
def rE():
    _ret = None
    _ret = 333.4
    return _ret

# Avogadros Number
def Avogadro():
    _ret = None
    _ret = 8.31441
    return _ret

# Gravitational constant
def GForce():
    _ret = None
    _ret = 9.81
    return _ret

# Kolar Mass of Air
def ML():
    _ret = None
    _ret = 28.96
    return _ret

# Molar Mass of Water
def MW():
    _ret = None
    _ret = 18.02
    return _ret

# 0 Degrees Celsius in K
def K0():
    _ret = None
    _ret = 273.15
    return _ret

# Pressure of water at the triple point in bar
def pwtr():
    _ret = None
    _ret = 0.00611657
    return _ret

# Temperature of Water at the triple point
def Twtr():
    _ret = None
    _ret = 273.16
    return _ret

# Temperature of water at the critical point
def Twcrit():
    _ret = None
    _ret = 647.096
    return _ret

def EnthalpyDry_p_rho(p_bar, rho):
    #
    _ret= (p_bar*100000*cpA())/(rho*RAir())-cpA() * K0() 
    return _ret

def Enthalpy_p_t_x(p_bar, t_Celsius, x_kgkg):
    _ret = None
    #
    # Calculates the enthlapy of humid air using
    # Airpressure
    # Temperature in Celsius
    # Absolute humidity in kg/kg
    #
    # It is calculated differently for saturation, wet area and ice area
    #
    pws = SaturationPressure_t(t_Celsius)
    #
    xs_l = HumidityRatioSat_t_p(t_Celsius, p_bar)
    #
    if  ( x_kgkg <= xs_l ) :
        #
        # ungesaettigte oder gerade gesaettigte feuchte Air
        #
        _ret = cpA() * t_Celsius + x_kgkg *  ( R0() + cpS() * t_Celsius )
        #
    else:
        if t_Celsius >= 0:
            #
            # Nebelgebiet, uebersaettigte feuchte Air
            #
            _ret = cpA() * t_Celsius + xs_l *  ( R0() + cpS() * t_Celsius )  +  ( x_kgkg - xs_l )  * cw() * t_Celsius
            #
        else:
            #
            # Eisgebiet
            #
            _ret = cpA() * t_Celsius + xs_l *  ( R0() + cpS() * t_Celsius )  +  ( x_kgkg - xs_l )  *  ( rE() - cE() * t_Celsius )
            #
    return _ret

def Enthalpy_p_t_phi(p_bar, t_Celsius, phi):
    _ret = None
    #
    # Calculates the enthlapy of humid air using
    # Airpressure
    # Temperature in Celsius
    # Relative humidity in 0 ... 0.5 = 50%
    #
    # It is calculated differently for saturation, wet area and ice area
    #
    pws = SaturationPressure_t(t_Celsius)
    #
    xs_l = HumidityRatioSat_t_p(t_Celsius, p_bar)
    #
    x_kgkg = HumidityRatio_p_t_phi(p_bar, t_Celsius, phi)
    #
    if  ( x_kgkg <= xs_l ) :
        #
        # ungesaettigte feuchte Air
        #
        _ret = cpA() * t_Celsius + x_kgkg *  ( R0() + cpS() * t_Celsius )
        #
    else:
        if t_Celsius >= 0:
            #
            # Nebelgebiet
            #
            _ret = cpA() * t_Celsius + xs_l *  ( R0() + cpS() * t_Celsius )  +  ( x_kgkg - xs_l )  * cw() * t_Celsius
            #
        else:
            #
            # Eisgebiet
            #
            _ret = cpA() * t_Celsius + xs_l *  ( R0() + cpS() * t_Celsius )  +  ( x_kgkg - xs_l )  *  ( rE() - cE() * t_Celsius )
            #
    #
    return _ret

def Enthalpy_t_x(t_Celsius, x_kgkg):
    _ret = None
    #
    # Calculates the enthlapy of humid air of 100% relative humidity using
    # 
    # Temperature in Celsius
    # Absolute humidity in kg/kg
    #
    _ret = cpA() * t_Celsius + x_kgkg *  ( R0() + cpS() * t_Celsius )
    #
    #
    return _ret

def Enthalpy_p_phi_x(p_bar, phi, x_kgkg):
    _ret = None
    #
    #
    pws = x_kgkg /  ( RA_RS() + x_kgkg )  * p_bar / phi
    #
    t_Celsius = SaturationTemperature_p(pws)
    #
    _ret = Enthalpy_p_t_x(p_bar, t_Celsius, x_kgkg)
    #
    return _ret

def HumidityRatioSat_t_p(Temperatur_C, Airpressure_bar):
    _ret = None
    #
    #
    # Calculates the humidity ratio in kg/kg for saturated air
    #
    pws = SaturationPressure_t(Temperatur_C)
    #
    # gibt den Druck in der Einheit bar zurueck
    #
    # aus Klimatechnik xs=0,622*pws/(p-pws)
    #
    _ret = RA_RS() * pws /  ( Airpressure_bar - pws )
    #
    return _ret

def HumidityRatio_p_t_phi(Airpressure_bar, Temperatur_C, RelativeHumidity):
    _ret = None
    # Calculates the humidity ratio x in kg/kg
    # Berechnet den Watergehalt x in kg/kg aus
    #
    # Airpressure         bar     airpressure
    # Temperatur        C       Temperature in Celsius
    # RelativeHumidity   0...    Relative humidity 1=100%   
    #
    pws = SaturationPressure_t(Temperatur_C)
    #
    _ret = RA_RS() * RelativeHumidity * pws /  ( Airpressure_bar - pws * RelativeHumidity )
    #
    return _ret

def HumidityRatio_p_pws_phi(Airpressure_bar, SaturationPressure_bar, RelativeHumidity):
    _ret = None
    # Calculates the humidity ratio x in kg/kg
    # Berechnet den Watergehalt x in kg/kg aus
    #
    # Airpressure         bar     airpressure
    # Saturation pressure  bar     Saturation pressure
    # RelativeHumidity   0...    Relative humidity 1=100% 
    #
    _ret = RA_RS() * RelativeHumidity * SaturationPressure_bar /  ( Airpressure_bar - SaturationPressure_bar * RelativeHumidity )
    #
    return _ret

def HumidityRatioSat_p_pws(p_bar, pws_bar):
    _ret = None
    # Calculates the humidity ratio x in kg/kg using airpressure and steam partial
    # pressure for saturated air 100% RH
    # returns the pressure in bar
    # gibt den Druck in der Einheit bar zurueck
    #
    # aus Klimatechnik xs=0,622*pws/(p-pws)
    #
    _ret = RA_RS() * pws_bar /  ( p_bar - pws_bar )
    #
    return _ret

def SaturationTemperature_p(Dampfdruck_bar):
    _ret = None
    # Returns Degrees Celsius
    # Calculates the saturation Temperature of water from -50 to 200 Degrees Celsius
    # using various models for different areas
    # -50 ... + 60  C Antoine Gleichung laut Baehr/Kabelac / Thermodynamik nach T aufgeloest
    # darueber bis 200 C Gleichungen von Bernd Glueck / Stoffwerte
    # bei Bereichsueberschreitung ist der Rueckgabewert -300
    #
    #
    p_bar = Dampfdruck_bar
    p_pa = Dampfdruck_bar * 100000
    #
    if p_bar >= 0.000039 and p_bar < pwtr():
        #
        # -50 bis Tripelpunkt
        #
        _ret = ( ( -22.5129 * Twtr() )  /  ( math.log(p_bar / pwtr()) - 22.5129 ) )  - K0()
        #
    elif p_bar >= pwtr() and p_bar < 0.19:
        #
        # Tripelpunkt bis + 60  C
        #
        _ret = ( ( -4102.99 )  /  ( math.log(p_bar / pwtr()) - 17.2799 ) )  - 237.431
        #
    elif p_bar >= 0.19 and p_bar < 1.0132:
        #
        # Berechnung nach Bernd Glueck, Stoffwerte, 1.5
        #
        g1 = - 63.16113
        g2 = 5.36859 * math.log(p_pa)
        g3 = 0.973587 * math.log(p_pa) ** 2
        g4 = - 0.0738636 * math.log(p_pa) ** 3
        g5 = 0.00481832 * math.log(p_pa) ** 4
        #
        _ret = g1 + g2 + g3 + g4 + g5
        #
    elif p_bar >= 1.0132 and p_bar < 15.6:
        #
        # Berechnung nach Bernd Glueck, Stoffwerte, 1.8
        #
        g1 = - 228.146
        g2 = 31.97037 * math.log(p_pa)
        g3 = 1.153295 * math.log(p_pa) ** 2
        g4 = - 0.27847109 * math.log(p_pa) ** 3
        g5 = 0.01319026 * math.log(p_pa) ** 4
        #
        _ret = g1 + g2 + g3 + g4 + g5
        #
    else:
        #
        _ret = - 300
        #
    #
    return _ret

def Temperature_h_x(h_kJkg, x_kgkg):
    _ret = None
    # Returns Degrees Celsius
    # Calculates the Temperature from given enthalpy and humidity ratio
    #
    _ret = ( h_kJkg - x_kgkg * R0() )  /  ( cpA() + x_kgkg * cpS() )
    #
    return _ret

def PartialPressureWater_p_x(p_bar, x_kgkg):
    _ret = None
    # Returns the partial pressure of steam in the air
    # Inputs Air pressure, Humidity ratio in kg/kg
    #
    _ret = p_bar * x_kgkg /  ( RA_RS() + x_kgkg )
    #
    return _ret

def PartialPressureAir_p_x(p_bar, x_kgkg):
    _ret = None
    #
    # Returns the partial pressure of air in the air
    # Inputs Air pressure, Humidity ratio in kg/kg
    #
    _ret = p_bar * RA_RS() /  ( RA_RS() + x_kgkg )
    #
    return _ret

def SaturationPressure_t(t_Celsius):
    _ret = None
    #
    # returns the saturation pressure of water
    # input Temperature in degrees celsius
    #
    # -50 ... + 60  C Antoine Gleichung laut Baehr/Kabelac / Thermodynamik
    # darueber bis 200 C Gleichungen von Bernd Glueck / Stoffwerte
    #
    T = t_Celsius + K0()
    #
    if t_Celsius >= -50 and t_Celsius < 0.01:
        #
        _ret = pwtr() * math.exp(22.5129 *  ( 1 - Twtr() / T ))
        #
    elif t_Celsius >= 0.01 and t_Celsius < 60:
        #
        _ret = pwtr() * math.exp(17.2799 -  ( 4102.99 /  ( t_Celsius + 237.431 ) ))
        #
    elif t_Celsius >= 60 and t_Celsius < 100:
        #
        # Berechnung nach Bernd Glueck, Stoffwerte, 1.4
        #
        g1 = - 0.000191275
        g2 = 0.07258 * t_Celsius
        g3 = - 0.0002939 * t_Celsius ** 2
        g4 = 0.0000009841 * t_Celsius ** 3
        g5 = - 0.00000000192 * t_Celsius ** 4
        #
        _ret = 611 * math.exp(g1 + g2 + g3 + g4 + g5) / 100000
        #
    elif t_Celsius >= 100 and t_Celsius <= 200:
        #
        # Berechnung nach Bernd Glueck, Stoffwerte, 1.7
        #
        g1 = 0.00006
        g2 = 0.0713274 * t_Celsius
        g3 = - 0.0002581631 * t_Celsius ** 2
        g4 = 0.0000006311955 * t_Celsius ** 3
        g5 = - 7.167112E-10 * t_Celsius ** 4
        #
        _ret = 611 * math.exp(g1 + g2 + g3 + g4 + g5) / 100000
        #
    else:
        # Ausserhalb des Gueltigkeitsbereiches -1
        _ret = - 1
        #
    #
    return _ret

def RelativeHumidity_p_x_pws(p_bar, x_kgkg, pws_bar):
    _ret = None
    #
    # Airpressure         bar             oder Pa
    # Watergehlat x in kg/kg
    # Saturation pressure  bar             oder Pa
    #
    # Airpressure und Saturation pressure muessen in der gleichen einheit uebergeben werden
    #
    # Formel (18)
    #
    _ret = ( p_bar / pws_bar )  * x_kgkg /  ( RA_RS() + x_kgkg )
    #
    return _ret

def RelativeHumidity_p_t_x(p_bar, t_Celsius, x_kgkg):
    _ret = None
    #
    #
    pws_bar = SaturationPressure_t(t_Celsius)
    #
    # Airpressure         bar             oder Pa
    # Watergehlat x in kg/kg
    # Saturation pressure  bar             oder Pa
    #
    # Airpressure und Saturation pressure muessen in der gleichen einheit uebergeben werden
    #
    # Formel (18)
    #
    _ret = ( p_bar / pws_bar )  * x_kgkg /  ( RA_RS() + x_kgkg )
    #
    return _ret

def RelativeHumidity_t_tDew(t_Celsius, tDew_Celsius):
    _ret = None
    #
    # Baehr S. 290
    #
    #
    pwstau = SaturationPressure_t(tDew_Celsius)
    pws = SaturationPressure_t(t_Celsius)
    #
    _ret = pwstau / pws
    #
    return _ret

def RelativeHumidity_p_t_twb(p_bar, t_Celsius, twb_Celsius):
    _ret = None
    #
    # Wir iterieren
    #
    #
    xs_fk = HumidityRatioSat_t_p(twb_Celsius, p_bar)
    #
    h1 = Enthalpy_p_t_x(p_bar, twb_Celsius, xs_fk)
    #
    # Das x fuer die gesuchte Relative Feuchte ist <= xs_fk
    #
    xs_gesucht = xs_fk
    #
    while 1:
        xs_gesucht = xs_gesucht - 0.001
        h2 = Enthalpy_p_t_x(p_bar, t_Celsius, xs_gesucht)
        if h2 < h1:
            break
    #
    while 1:
        xs_gesucht = xs_gesucht + 0.00001
        h2 = Enthalpy_p_t_x(p_bar, t_Celsius, xs_gesucht)
        if h2 >= h1:
            break
    #
    _ret = RelativeHumidity_p_t_x(p_bar, t_Celsius, xs_gesucht)
    #
    return _ret

def WetBulbTemp_p_t_phi(p_bar, t_Celsius, phi):
    _ret = None
    #
    # Berechnet iterativ die Feuchtkugeltemperatur
    #
    #
    # Airpressure         bar
    # Temperatur         C
    # RelativeHumidity   0...1
    #
    # 1. Die enthalpie des Zustandspunktes wird berechnet
    #
    h1 = Enthalpy_p_t_phi(p_bar, t_Celsius, phi)
    #
    # 2.    Die Saettigungstemperatur mit der gleichen enthalpie wird gesucht
    #       ... so als wuerden wir der Isenthalpen auf einem hx-Diagramm
    #       zur Saettigungslinie hin folgen
    #
    twb = - 50
    #
    while 1:
        twb = twb + 1
        hfk = Enthalpy_p_t_phi(p_bar, twb, 1)
        if h1 < hfk:
            break
    #
    while 1:
        twb = twb - 0.001
        hfk = Enthalpy_p_t_phi(p_bar, twb, 1)
        if h1 >= hfk:
            break
    #
    #
    _ret = twb
    #
    #
    return _ret

def Density_p_t_pw(Airpressure_bar, Temperatur_C, Partialdampfdruck_bar):
    _ret = None
    #
    # Airpressure             bar
    # Partialdampfdruck     bar
    # Temperatur             C
    #
    # Direkt nach Baehr mit Verwendung der definierten Konstanten
    #
    # Zur Uebereinstimmung mit BFS Formel 19 muessen RL und RD angeglichen werden
    #
    #
    Kelvintemperatur = Temperatur_C + K0()
    Partialdampfdruck = Partialdampfdruck_bar * 100000
    Airpressure = Airpressure_bar * 100000
    #
    _ret = Airpressure /  ( Kelvintemperatur * RAir() )  *  ( 1 -  ( ( Partialdampfdruck / Airpressure )  *  ( 1 -  ( RAir() / RSteam() ) ) ) )
    #
    return _ret

def Density_p_t_x(p_bar, t_Celsius, x_kgkg):
    _ret = None
    #
    #
    Pw_bar = PartialPressureWater_p_x(p_bar, x_kgkg)
    #
    _ret = Density_p_t_pw(p_bar, t_Celsius, Pw_bar)
    #
    return _ret

def Density_p_t_phi(p_bar, t_Celsius, phi):
    _ret = None
    #
    #
    x_kgkg = HumidityRatio_p_t_phi(p_bar, t_Celsius, phi)
    #
    _ret = Density_p_t_x(p_bar, t_Celsius, x_kgkg)
    #
    return _ret

def DensitySat_p_t(p_bar,t_Celsius):
    #
    pws=SaturationPressure_t(t_Celsius)*100000
    #
    pressure=p_bar*100000
    T_Kelvin=t_Celsius+K0()
    #
    _ret=(((pressure-pws)/RAir())+(pws/RSteam()))/T_Kelvin
    return _ret

def TemperatureSat_p_rho(p_bar,rho):
    #
    # Returns Celsius Temperature for 100%, at given pressure and density
    #
    _ret=brentq( lambda t: DensitySat_p_t(p_bar,t)-rho,-50,200,xtol=0.001)
    return _ret
    
def EnthalpySat_p_rho(p_bar,rho):
    #
    t_Celsius=TemperatureSat_p_rho(p_bar,rho)
    #
    _ret=Enthalpy_p_t_phi(p_bar, t_Celsius, 1)
    #
    return _ret
    
def HumidityRatioSat_p_rho(p_bar,rho):
    #
    t_Celsius=TemperatureSat_p_rho(p_bar,rho)
    #
    _ret=HumidityRatioSat_t_p(t_Celsius, p_bar)
    return _ret    
    
def Reynolds_d_nue_w(di_m, nue_m2s, w_ms):
    _ret = None
    #
    _ret = w_ms * di_m / nue_m2s
    #
    return _ret

def Reynolds_di_eta_rho_w(di_m, eta_Pas, rho_kgm3, w_ms):
    _ret = None
    #
    _ret = rho_kgm3 * w_ms * di_m / eta_Pas
    #
    return _ret

def Prandtl_t(t_Celsius):
    _ret = None
    #
    # Nach Bernd Glueck, Stoffwerte (2.45)
    #
    #
    g1 = 0.71789
    g2 = - 0.0001675739 * t_Celsius
    g3 = 0.0000006514142 * t_Celsius ** 2
    g4 = - 6.687762E-10 * t_Celsius ** 3
    #
    _ret = g1 + g2 + g3 + g4
    return _ret

def ViskosityKin_t(t_Celsius):
    _ret = None
    #
    # Nach Bernd Glueck, Stoffwerte (2.41)
    # in Nue m^2/s fuer 1 bar
    #
    # RM 9.7.12 Die abweichungen gegenueber Refprop sind mit
    # maximal 0,412% bei -50 C bis -20 C ab dort bis + 200 C MPa maximaler Fehler
    # gegenueber Refprop 0,125%
    #
    #
    g1 = 0.0000135198
    g2 = 0.00000008930841 * t_Celsius
    g3 = 1.094808E-10 * t_Celsius ** 2
    g4 = - 3.659345E-14 * t_Celsius ** 3
    #
    _ret = g1 + g2 + g3 + g4
    #
    return _ret

def ViskosityDyn_t(t_Celsius):
    _ret = None
    #
    # Nach Bernd Glueck, Stoffwerte (2.40)
    # Eta in Pas bei einem Druck von 1 bar
    # RM 9.7.12 Die Abweichungen gegenueber Refprop sind mit
    # maximal 0,22% bei 0.1 MPa und -50 .. + 200 C MPa gering
    #
    #
    g1 = 0.0000172436
    g2 = 0.0000000504587 * t_Celsius
    g3 = - 3.923361E-11 * t_Celsius ** 2
    g4 = 4.046118E-14 * t_Celsius ** 3
    #
    _ret = g1 + g2 + g3 + g4
    #
    return _ret

def Waermeleitfaehigkeit_t(t_Celsius):
    _ret = None
    #
    # in W/(mK)
    #
    #
    g1 = 0.024178
    g2 = 0.00007634878 * t_Celsius
    g3 = - 0.00000004663859 * t_Celsius ** 2
    g4 = 4.612639E-11 * t_Celsius ** 3
    #
    _ret = g1 + g2 + g3 + g4
    #
    return _ret

def FrictionFactor(Reynolds, Rauhigkeit_m, Diameter_m):
    _ret = None
    #
    #   (c) Reiner Johannes Mayers 2010
    #       Version 0.1 Januar 2010
    #
    #   Quellenangaben:
    #   Stroemung und Druckverlust, Kamprath
    #   Pohlmann, Handbuch der Kaeltetechnik
    #   Wikipedia engl. Druckverlust, FrictionFactor (Serghides solution)
    #
    #
    rrk = Rauhigkeit_m / Diameter_m
    rrk37 = rrk / 3.7
    #
    #
    if Reynolds >= 0 and Reynolds <= 2320:
        #
        # FrictionFactor nach Hagen Poisseuille fuer laminare
        # Stroemung unter Vernachlaessigung der RohrRauhigkeit_m
        #
        _ret = 64 / Reynolds
        #
        #ElseIf Reynolds > 2320# And Reynolds <= 2500# Then
        #
        # FrictionFactor nach Blasius fuer turbulente
        # Stroemung unter Vernachlaessigung der RohrRauhigkeit_m
        #
        #FrictionFactor = 0.3164 / (Reynolds ^ (0.25))
        #
    elif Reynolds > 2320:
        #
        # FrictionFactor nach Serghides Naeherungsformel fuer die Colebrook White
        # Formel unter Beruecksichtigung der RohrRauhigkeit_m
        #
        sga = - 2 * math.log10(rrk37 +  ( 12 / Reynolds ))
        sgb = - 2 * math.log10(rrk37 +  ( ( 2.51 * sga )  / Reynolds ))
        sgc = - 2 * math.log10(rrk37 +  ( ( 2.51 * sgb )  / Reynolds ))
        _ret = 1 /  ( sga -  ( ( ( sgb - sga )  ** 2 )  /  ( sgc -  ( 2 * sgb )  + sga ) ) )  ** 2
    return _ret

def FrictionFactorDC(Reynolds, Rauhigkeit_m, Diameter_m):
    _ret = None
    #
    #   (c) Reiner Johannes Mayers 2010
    #       Version 0.1 Maerz 2010
    #
    #   Quellenangaben:
    #   Didier Clamond, Efficient Resolution of the Colebrook Equation
    #   American Chemical Society, Published on Web 02/13/2009
    #   Die Version wurde aus dem Fortran Code der obigen Publikation erstellt
    #
    #   Eine vergleichbare Praezision wie in Fortran ist mit VB ohne definition eigener
    #   Rechenoperationen nicht moeglich VB kuerzt die gegebenen Konstanten um 4 Stellen
    #   C1 =                 1.151292546497022842
    #   X1 = K * Reynolds *  0.123968186335417556
    #   X2 = Log(Reynolds) - 0.779397488455682028
    #
    #
    if Reynolds >= 0 and Reynolds <= 2300:
        #
        # FrictionFactor nach Hagen Poisseuille fuer laminare
        # Stroemung unter Vernachlaessigung der RohrRauhigkeit_m
        #
        _ret = 64 / Reynolds
        #
        #ElseIf Reynolds > 2320# And Reynolds <= 2500# Then
        #
        # FrictionFactor nach Blasius fuer turbulente
        # Stroemung unter Vernachlaessigung der RohrRauhigkeit_m
        #
        #FrictionFactor = 0.3164 / (Reynolds ^ (0.25))
        #
    elif Reynolds > 2300: 
        #
        K = Rauhigkeit_m / Diameter_m
        T = 1 / 3
        #
        # Initialisierung
        #
        x1 = K * Reynolds * 0.123968186335418
        x2 = math.log(Reynolds) - 0.779397488455682
        C1 = 1.15129254649702
        #
        # Anfangswert Schaetzung
        #
        F = x2 - 0.2
        #
        # Erste Iteration
        #
        E = ( math.log(x1 + F) - 0.2 )  /  ( 1 + x1 + F )
        F = F -  ( ( 1 + x1 + F + 0.5 * E )  * E *  ( x1 + F ) )  /  ( 1 + x1 + F +  ( E *  ( 1 + E * T ) ) )
        #
        # Zweite Iteration wenn noetig
        #
        if  ( ( x1 + x2 )  < 5.7 ) :
            E = ( math.log(x1 + F) + F - x2 )  /  ( 1 + x1 + F )
            F = F -  ( ( 1 + x1 + F + 0.5 * E )  * E *  ( x1 + F ) )  /  ( 1 + x1 + F +  ( E *  ( 1 + E * T ) ) )
            #
            # Endgueltige Berechnung
            #
            F = C1 / F
            _ret = F * F
            #
    return _ret

def DeltaP(Lambda, L_m, d_m, rho_kgm3, w_ms):
    _ret = None
    #
    # Bernoulli
    #
    _ret = ( Lambda * L_m * rho_kgm3 * w_ms ** 2 )  /  ( d_m * 2 )
    #
    return _ret

def DewPointTemperature(pwd_bar):
    _ret = None
    _ret = SaturationTemperature_p(pwd_bar)
    return _ret

def HeatCapacityAir_t(t_Celsius):
    _ret = None
    #
    # Calculates cp of dry air for atmospheric pressures
    #
    # Berechnet die spezifische Waermekapazitaet der trockenen Air
    # ... eigentlich fuer genau 1 bar
    # Beim Vergleich mit Refprop, Air mit 800 mBar und 1100 mBar bleiben
    # die Abweichungen unter 0,2%, bei 1bar sogar unter 0,1%
    #
    #
    if t_Celsius < 0:
        # Im Minusbereich wird eine Formel aus dem Handbuch der Software AHH Verwendet
        # Gegenueber den Formeln von Glueck zeigt sie geringere Abweichungen zu Refprop
        #
        A = 5.240699E-13
        B = - 2.979734E-10
        C = 0.000000172963
        D = 0.0000101051
        E = 1.006156
        #
        T = t_Celsius
        #
        _ret = A * T ** 4 + B * T ** 3 + C * T ** 2 + D * T + E
    else:
        # Im Plusbereich passen die Formeln von Glueck besser
        # im Vergleich zu Refprop
        #
        g1 = 1.0065
        g2 = 0.000005309587 * t_Celsius
        g3 = 0.0000004758596 * t_Celsius ** 2
        g4 = - 1.136145E-10 * t_Celsius ** 3
        #
        _ret = g1 + g2 + g3 + g4
        #
    #
    return _ret

def HeatCapacitySteam_t(t_Celsius):
    _ret = None
    #
    # Reiner Mayers 11.7.2012
    #
    # Im Minusbereich stimmen die Werte nach Baehr ganz gut mit den mir bekannten
    # Tabellenwerten ueberein, im Plusbereich kommen steigend schnell grosse
    # Abweichungen gegenueber Refprop. Daher wird ab 5 C auf die Berechnungsgleichung
    # nach Glueck umgestellt Stoffwerte 1.45.
    # bis 50 C ergeben sich so Abweichungen kleiner 1% gegenueber Refprop (IAWPS97)
    # bis 200 C steigt die Abweichung auf bis zu 3,75%
    #
    # Berechnung der Waermekapazitaet von Dampf nach
    # Baehr,H.D.; Diederichsen, Chr.: "Berechnungsgleichungen fuer Enthalpie und
    # Entropie der Komponenten von Air und Verbrennungsgasen 1988
    #
    #
    #
    pk[1] = - 0.0020169
    pk[2] = 0.05324806
    pk[3] = - 0.59025452
    pk[4] = 3.57748654
    pk[5] = - 12.97778482
    pk[6] = 33.32045416
    pk[7] = - 42.64731893
    pk[8] = 41.39089229
    pk[9] = - 23.96053283
    pk[10] = 8.22692113
    pk[11] = - 1.55193837
    pk[12] = 0.12417097
    #
    if t_Celsius <= 5:
        # nach Baehr
        #
        T = ( t_Celsius + K0() )  / 1000
        #
        cpD = 0
        for Zaehler in vbForRange(1, 12):
            cpD = cpD + pk(Zaehler) * T **  ( Zaehler - 6 )
        cpD = cpD * RSteam() / 1000
        #
        _ret = cpD
        #
    else:
        # nach Glueck
        #
        g1 = 1.854283
        g2 = 0.00112674 * t_Celsius
        g3 = - 0.000006939165 * t_Celsius ** 2
        g4 = 0.0000001344783 * t_Celsius ** 3
        #
        _ret = g1 + g2 + g3 + g4
    #
    return _ret

def MerkelsNumber(t_WaterEin_Celsius, t_WaterAus_Celsius, t_AirEin_Celsius, t_FeuchtkugelEin_Celsius, Airzahl, Airpressure_bar, N=1000):
    _ret = None
    #
    # Berechnung grob angelehnt an VDI Energietechnische Arbeitsmappe 14. Aufl. 1995
    # Es werden die vorhandenen Rechenroutinen fuer Saettigungsdruecke verwendet
    #
    #
    #
    _ret = 0
    J = 0
    Mer = 0
    #
    # 1. RelativeHumidity_p_t_twb, recht genau ueber eine iteration
    # ... darueber dann x und h
    #
    phi = RelativeHumidity_p_t_twb(Airpressure_bar, t_AirEin_Celsius, t_FeuchtkugelEin_Celsius)
    xL = HumidityRatio_p_t_phi(Airpressure_bar, t_AirEin_Celsius, phi)
    hL = Enthalpy_t_x(t_AirEin_Celsius, xL)
    A = ( t_WaterEin_Celsius - t_WaterAus_Celsius )  / N
    #
    while 1:
        #
        tW = t_WaterAus_Celsius + A *  ( 0.5 + J )
        #
        # Nur die Saettigungszustaende
        #
        xS = HumidityRatioSat_t_p(tW, Airpressure_bar)
        hLS = Enthalpy_t_x(tW, xS)
        hm = hL + cw() * A / Airzahl *  ( 0.5 + J )
        K = cw() * A /  ( hLS - hm )
        Mer = Mer + K
        J = J + 1
        if ( J == N - 1 ):
            break
    #
    _ret = Mer
    #
    return _ret


def BarometricPressure(hoehe_H, Airfeuchte_phi, t_jahresmittel):
    _ret = None
    #
    #
    H = hoehe_H
    T = t_jahresmittel + Twtr()
    phi = Airfeuchte_phi
    #
    # Zuerst ohne Airfeuchte
    #
    p1 = 1.01325 * 1 / math.exp(( ( ML() * GForce() * H )  /  ( Avogadro() * 1000 * T ) ))
    pd = SaturationPressure_t(t_jahresmittel)
    x = ( MW() / ML() )  *  ( phi * pd )  /  ( p1 - phi * pd )
    #
    _ret = 1.01325 * 1 / math.exp(( ( ML() * GForce() * H )  /  ( Avogadro() * 1000 * T ) )  *  ( ( 1 + x )  /  ( 1 + x * ML() / MW() ) ))
    #
    return _ret

def MassFlowDry(Vfl_m3h, rho, x_kgkg):
    _ret = None
    #
    # Massenstrom der trockenen Air
    #
    _ret = Vfl_m3h * rho /  ( 1 + x_kgkg )
    #
    return _ret

def VolumeFlowHumid(mstromtr_kgh, rho, x_kgkg):
    _ret = None
    #
    # Volumenstrom der feuchten Air in m^3/h
    #
    _ret = mstromtr_kgh *  ( 1 + x_kgkg )  / rho
    #
    return _ret

def VolumeFlowm3h_A_w(A_m2, w_ms):
    _ret = None
    #
    _ret = A_m2 * w_ms * 3600
    #
    return _ret

def Velocity_A_VFlow(A_m2, vstrom_m3h):
    _ret = None
    #
    #
    vstrom_m3s = vstrom_m3h / 3600
    #
    _ret = vstrom_m3s / A_m2
    #
    return _ret

def Area_VFlow_w(vstrom_m3h, w_ms):
    _ret = None
    #
    #
    vstrom_m3s = vstrom_m3h / 3600
    #
    _ret = vstrom_m3s / w_ms
    #
    return _ret

def AreaCircle(Diameter_m):
    _ret = None
    #
    _ret = Diameter_m ** 2 * math.pi / 4
    #
    return _ret

def AreaSquare(b_m, h_m):
    _ret = None
    #
    # Querschnittsflaeche eines Rechteckigen Airkanals
    # mit hoehe h und Breite b in m
    #
    _ret = b_m * h_m
    #
    return _ret

def HydraulicDiameter(b_m, h_m):
    _ret = None
    #
    # Berechnet den hydraulischen Diameter eines Rechteckigen Airkanals
    # mit hoehe h und Breite b in m
    #
    _ret = 4 * b_m * h_m /  ( 2 * b_m + 2 * h_m )
    #
    return _ret

def Diameter(A_m2):
    _ret = None
    #
    _ret = math.sqrt(A_m2 * 4 / math.pi())
    #
    return _ret

def RelativeHumidityHeat(p_bar, t_Celsius_Ein, T_Celsius_Aus, phi_Ein):
    _ret = None
    #
    # Aufheizen der Air bei gleichem Airpressure, x= const.
    #
    #
    x = HumidityRatio_p_t_phi(p_bar, t_Celsius_Ein, phi_Ein)
    #
    _ret = RelativeHumidity_p_t_x(p_bar, T_Celsius_Aus, x)
    #
    return _ret

def MixX(x1, m1, x2, m2):
    _ret = None
    _ret = ( x1 * m1 + x2 * m2 )  /  ( m1 + m2 )
    return _ret

def MixH(h1, m1, h2, m2):
    _ret = None
    _ret = ( h1 * m1 + h2 * m2 )  /  ( m1 + m2 )
    return _ret

def MixT(t1, m1, t2, m2):
    _ret = None
    _ret = ( t1 * m1 + t2 * m2 )  /  ( m1 + m2 )
    return _ret
