import scipy.special as special
import math
import matplotlib.pyplot as plt

cte = ( 2 / math.pi ) ** 2
K = 10                                   # unit = md
h = 150                                  # unit = ft
q = 1000                                 # unit = STB / D
B = 1.475                                # unit = RB / STB
miu = .72                                # unit = cp
ct = 1.5 * 10 ** ( - 5 )                 # unit = 1 / psi
rw = .5                                  # unit = ft
phi = .23

s_List = [ 0 , 5 , 10 , 20 ]
CD_List = [ 0 , 10**2 , 10**3 , 10**4 , 10**5]
tD_List = [ 10 ** 2 , 2 * 10 ** 2 , 3 * 10 ** 2 , 5 * 10 ** 2 , 7 * 10 ** 2 , 10 ** 3 , 2 * 10 ** 3 , 3 * 10 ** 3 , 5 * 10 ** 3 , 7 * 10 ** 3 , 10 ** 4 , 2 * 10 ** 4 , 3 * 10 ** 4 , 5 * 10 ** 4 , 7 * 10 ** 4 , 10 ** 5 , 2 * 10 ** 5 , 3 * 10 ** 5 , 5 * 10 ** 5 , 7 * 10 ** 5 , 10 ** 6 , 2 * 10 ** 6 , 3 * 10 ** 6 , 5 * 10 ** 6 , 7 * 10 ** 6 , 10 ** 7 , 2 * 10 ** 7 , 3 * 10 ** 7 , 5 * 10 ** 7 , 7 * 10 ** 7 , 10 ** 8 ]
C_List = []
T_List = []

P_factor = ( 141.2 * q * miu * B  ) / ( K * h )

for CD in CD_List :
    C = ( CD * phi * ct * h * rw ** 2 ) / .8936
    C_List.append( C )

for tD in tD_List :
    t = ( tD * phi * miu * ct * rw ** 2 ) / ( .0002637 * K )
    T_List.append( t )


def j0 ( x ) :
    return special.jv( 0 , x )

def j1 ( x ) :
    return special.jv( 1 , x )

def y0 ( x ) :
    return special.yv( 0 , x )

def y1 ( x ) :
    return special.yv( 1 , x )

def f ( x , tD , CD , s ) :
    return  cte * ( ( 1 - math.exp ( - x ** 2 * tD ) ) / ( x ** 3 * ( ( x * CD * j0 ( x ) - ( 1 - CD * s * x ** 2 ) * j1 ( x ) ) ** 2 + ( x * CD * y0( x ) - ( 1 - CD * s * x ** 2 ) * y1 ( x ) ) ** 2 ) ) )

for i in range ( len ( CD_List ) ) :
    CD = CD_List [ int ( i ) ]
    C = C_List [ int ( i ) ]

    for s in s_List :
        PwD = []
        delta_P = []

        for tD in tD_List :

            def simpson ( f , a , b , m ) :
                k = 0
                h = ( b - a ) / m

                for i in range ( 1 , m ) :
                    if i % 3 == 0 :
                        k += 2 * f ( a + i * h , tD , CD , s )
                    else :
                        k += 3 * f ( a + i * h , tD , CD , s )

                k += f ( b , tD , CD , s )
                PD = ( 3 * h * k ) / 8
                DELTA_P = P_factor * PD
                PwD.append( PD )
                delta_P.append( DELTA_P )

            simpson ( f , 0 , .2 , 10000 )

        plt.figure( 'The Ramey type curve __ Dimensionless' , figsize = ( 19 , 9 ) )
        plt.plot( tD_List , PwD , label = 'CD = {} , s = {}'.format ( CD , s ) )
        plt.xlabel( 't_D' )
        plt.ylabel( 'P_D' )
        plt.xscale( 'log' )
        plt.yscale( 'log' )
        plt.legend()
        plt.tight_layout()

        plt.figure( 'The Ramey type curve' , figsize = ( 19 , 9 ) )
        plt.plot( T_List , delta_P , label = 'C = {} , s = {}'.format( C , s ) )
        plt.xlabel(' t (hours) ')
        plt.ylabel(' delta_P (psi) ')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
        plt.tight_layout()
        
plt.show()
