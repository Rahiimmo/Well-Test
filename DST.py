import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl

#inputs
P = [12.10762332,30.2690583,84.75336323,181.6143498,272.4215247,339.0134529,447.9820628,611.4349776,732.5112108,871.7488789,1023.09417,1132.06278,1247.085202,1434.753363,1537.668161,1598.206278,1695.067265,1779.820628,1888.789238,1961.434978,2003.811659,2009.865471,2003.811659,2028.026906,1997.757848,1755.605381,1580.044843,478.2511211,508.5201794,544.8430493,581.1659193,611.4349776,641.7040359,678.0269058,708.2959641,780.941704,889.9103139,980.7174888,1029.147982,1089.686099,1132.06278,1162.331839,1174.439462,1192.600897,1204.70852,1210.762332,1210.762332,1216.816143,1216.816143,1144.170404,944.3946188,829.3721973,732.5112108,762.7802691,786.9955157,817.264574,841.4798206,871.7488789,902.0179372,938.3408072,968.6098655,998.8789238,1023.09417,1053.363229,1065.470852,1089.686099,1113.901345,1126.008969,1144.170404,1192.600897,1404.484305,1695.067265,1840.358744,2003.811659,2221.748879,2488.116592,2494.170404,2518.38565,2536.547085,2572.869955,2385.201794,2034.080717,2021.973094,1973.542601,1785.874439,1555.829596,1374.215247,1071.524664,823.3183857,556.9506726,339.0134529,175.5605381,60.53811659,24.21524664,24.21524664]
t = [0,0.86381323,1.26692607,1.785214008,2.764202335,3.858365759,4.203891051,4.66459144,4.952529183,5.298054475,5.643579767,5.93151751,6.219455253,6.737743191,6.910505837,7.083268482,7.313618677,7.486381323,7.889494163,8.062256809,8.638132296,9.501945525,12.72684825,13.41789883,14.10894942,14.22412451,14.23547863,14.56964981,15.03035019,15.66381323,16.41245136,16.98832685,17.56420233,18.42801556,19.1766537,19.46459144,19.63735409,19.92529183,20.21322957,20.90428016,21.65291829,23.03501946,24.24435798,26.49027237,27.64202335,29.25447471,30.80933852,31.96108949,32.99766537,33.11284047,33.21547896,33.41335977,33.51660449,33.68871595,34.14941634,34.95564202,35.81945525,36.74085603,38.06536965,39.44747082,41.34785992,43.19066148,44.80311284,46.5307393,48.37354086,50.27392996,52.98054475,55.28404669,57.1844358,57.2264907,57.29961089,57.35719844,57.58754864,57.70272374,58.16342412,58.45136187,60.29416342,61.10038911,61.50350195,62.02178988,62.07937743,62.19455253,63.05836576,63.97976654,64.78599222,65.41945525,65.82256809,66.85914397,67.492607,68.2412451,468.81712062,69.7385214,70.31439689,71.00544747,71.69649805]
q = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.96528,27.309,23.2986,20.816,19.0972,18.33333,15.85069444,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.71875,11.7468,16.5346,21.5799,20.816,16.9965,13.3681,12.2222,9.93056,9.54861,8.40278,6.49306,6.11111,4.58333,4.39236,3.81944,3.24653,3.05556,2.29167,1.90972,0.84567,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
# Units:        P >> psi       ,    t >> hr       ,    q >> STBD
h = 50                                  # unit = ft
B = 1.475                               # unit = RB/STB
miu = 1.72                              # unit = cp
ct = 1.5 * 10**(-5)                     # unit = 1/psi
rw = 0.5                                # unit = ft
phi = 0.23

draw_down_1_t = []                      # time data for the first drawdown 
draw_down_1_q = []                      # flow rate data for the for the first drawdown
draw_down_2_t = []                      # time data for the second drawdown
draw_down_2_q = []                      # flow rate data for the for the second drawdown
build_up_1_P = []                       # Pressure data for the first build-up 
build_up_1_t = []                       # time data for the first build-up
build_up_2_P = []                       # Pressure data for the second build-up
build_up_2_t = []                       # time data for the second build-up

# Identifying the first and second drawdowns 
for i in range( 0 , len ( q ) ) :

    if q [ i ] != 0 :
        if len ( draw_down_1_t ) == 0 :
            draw_down_1_t.append( t [ i ] )
            draw_down_1_q.append( q [ i ] )
            start_dd1_t = t [ i - 1 ]                 # t0 for the first drawdown

        else :
            if q [ i - 1 ] != 0 :
                if len ( draw_down_2_t ) == 0 :
                    draw_down_1_t.append( t [ i ] )
                    draw_down_1_q.append( q [ i ] )

                elif len ( draw_down_2_t ) != 0 :
                    draw_down_2_t.append( t [ i ] )
                    draw_down_2_q.append( q [ i ] )

            else :
                draw_down_2_t.append( t [ i ] )
                draw_down_2_q.append( q [ i ] )
                start_dd2_t = t [ i - 1 ]           # t0 for the second drawdown

# Identifying the first and second build-ups 
j = True
for i in range ( 0, len ( q ) ) :

    if j == True :
        if draw_down_1_t [ len ( draw_down_1_t ) - 1 ] < t [ i ] < draw_down_2_t [ 0 ] :
            build_up_1_t.append( t [ i ] )
            build_up_1_P.append( P [ i ] )

        else:
            if draw_down_2_t [ len ( draw_down_2_t ) - 1 ] < t [ i ] :
                if P [ i ] - P [ i - 1 ] > 10 :
                    build_up_2_t.append( t [ i ] )
                    build_up_2_P.append( P [ i ] )

                else :
                    build_up_2_t.append( t [ i ] )
                    build_up_2_P.append( P [ i ] )
                    j = False

# Calculation of t0 and its related pressure for the first and second build-ups  
for i in range ( 1, len ( P) ) :

    if P [ i ] == build_up_1_P [ 0 ] :
        start_bu1_P = P [ i - 1 ]
        start_bu1_t = t [ i - 1 ]

    if P [ i ] == build_up_2_P [ 0 ] :
        start_bu2_P = P [ i - 1 ]
        start_bu2_t = t [ i - 1 ]

#update the lists
dd_1_t = [0]             # t0 for odeh , selig formula
for i in draw_down_1_t :
    f = i - start_dd1_t
    dd_1_t.append( f )

dd_1_q = [0] + draw_down_1_q

dd_2_t = [0]
for i in draw_down_2_t :
    g = i - start_dd2_t
    dd_2_t.append( g )

dd_2_q = [0] + draw_down_2_q

bu_1_t = [0]
for i in build_up_1_t :
    z = i - start_bu1_t
    bu_1_t.append( z )

bu_1_P = [ start_bu1_P ] + build_up_1_P

bu_2_t = [0]
for i in build_up_2_t :
    x = i - start_bu2_t
    bu_2_t.append( x )

bu_2_P = [ start_bu2_P ] + build_up_2_P

#odeh and selig
a1 = 0                 # Soorat-e formule
a2 = 0                 # Makhraj-e formule
for  i in range ( 1 , len ( dd_1_t ) ) :
    a1 += dd_1_q [ i ] * ( ( dd_1_t [ i ])**2 - ( dd_1_t [ i - 1] )**2 )
    a2 += dd_1_q [ i ] * ( dd_1_t [ i ] - dd_1_t [ i - 1] )

tp1 = 2 * ( dd_1_t [ len ( dd_1_t ) - 1 ] - a1 / ( a2 * 2 ) )
q1 = a2 / tp1

a3 = 0                  # The numerator
a4 = 0                  # The denominator
for  i in range ( 1 , len ( dd_2_t ) ) :
    a3 += dd_2_q [ i ] * ( ( dd_2_t [ i ] )**2 - ( dd_2_t [ i - 1 ] )**2 )
    a4 += dd_2_q [ i ] * ( dd_2_t [ i ] - dd_2_t [ i - 1] )

tp2 = 2 * ( dd_2_t [ len ( dd_2_t ) - 1 ] - a3 / ( a4 * 2 ) )
q2 = a4 / tp2

#hornor time
hornor1 = []
for i in range ( 1 , len ( bu_1_t ) ) :
    j = bu_1_t [ i ]
    h1 = ( j + tp1 ) / j
    hornor1.append( h1 )

hornor2 = []
for i in range ( 1 , len ( bu_2_t ) ) :
    j = bu_2_t [ i ]
    h2 = ( j + tp2 ) / j
    hornor2.append( h2 )

plt.figure( 'Build up 1' )
plt.scatter( hornor1 , build_up_1_P )
plt.xlabel( 'hornor time' )
plt.ylabel( 'P_ws' )
plt.xscale( 'log' )

plt.figure( 'Build up 2' )
plt.scatter( hornor2 , build_up_2_P )
plt.xlabel( 'hornor time' )
plt.ylabel( 'P_ws' )
plt.xscale( 'log' )


# trendline 1 calculation
z = np.polyfit( np.log10 ( hornor1 ) , build_up_1_P , 1 )
p = np.poly1d( z )
pl.plot( x , p( x ) , "r--" )
# the line equation :
print( "P_ws=%.6f*Log(x)+(%.6f)"%( z [ 0 ] , z [ 1 ] ) )

# trendline 2 calculation
w = np.polyfit( np.log10 ( hornor2 ) , build_up_2_P , 1 )
p = np.poly1d( w )
pl.plot( x , p(x) , "r--" )
# the line equation :
print( "P_ws=%.6f*Log(x)+(%.6f)"%( w [ 0 ] , w [ 1 ] ) )

m1 = z [ 0 ]                                 # slope of the first plot
m2 = w [ 0 ]                                 # slope of the second plot

k1 = 162.6 * q1 * B * miu / ( abs ( m1 ) * h )
k2 = 162.6 * q1 * B * miu / ( abs ( m2 ) * h )
P_1hr_1 = m1 * np.log10( tp1 + 1 ) + z [ 1 ]
P_1hr_2 = m2 * np.log10( tp2 + 1 ) + w [ 1 ]
S1 = 1.151 * ( ( P_1hr_1 - start_bu1_P ) / abs ( m1 ) - np.log10( k1 / ( phi * miu * ct * ( rw**2 ) ) ) + 3.23 + np.log10( ( tp1 + 1 ) / tp1 ) )
S2 = 1.151 * ( ( P_1hr_2 - start_bu2_P ) / abs ( m2 ) - np.log10( k2 / ( phi * miu * ct * ( rw**2 ) ) ) + 3.23 + np.log10( ( tp2 + 1 ) / tp2 ) )
P_star_1 = z [ 1 ]
P_star_2 = w [ 1 ]

print( 'time_drawdown1 = ' , draw_down_1_t )
print( 'flow_rate_drawdown1 = ' , draw_down_1_q )
print( 'time_drawdown2 = ' , draw_down_2_t )
print( 'flow_rate_drawdown2 = ' , draw_down_2_q )
print( 'time_build_up1 =  ' , build_up_1_t )
print( 'Pressure_build_up1 = ' , build_up_1_P )
print( 'time_build_up2 =  ' , build_up_2_t )
print( 'Pressure_build_up2 = ' , build_up_2_P )

print( 'k_1 =' , k1 , 'md' )
print( 'k_2 =' , k2 , 'md' )
print( 'S_1 =' , S1 )
print( 'S_2 =' , S2 )
print( 'P_star_1 =' , P_star_1 , 'psi' )
print( 'P_star_2 =' , P_star_2 , 'psi' )

plt.tight_layout()
plt.show()
