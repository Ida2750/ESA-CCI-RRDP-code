def SnowDepth(lat, lon, month):

    from numpy import sin, cos, pi

    Y = (90-lat) * sin(lon*(pi/180))
    X = (90-lat) * cos(lon*(pi/180))


    #Set the coefficients

    if month == 1: #JANUARY
        H_0 = 28.01
        A = 0.1270
        B = -1.1833
        C = -0.1164
        D = -0.0051
        E = 0.0243
        epsilon = 7.6
        F = -0.06
        sigma_f = 0.07
        IAV = 4.6

    elif month == 2: #FEBRUARY
        H_0 = 30.28
        A = 0.1056
        B = -0.5908
        C = -0.0263
        D = -0.0049
        E = 0.0044
        epsilon = 7.9
        F = -0.06
        sigma_f = 0.08
        IAV = 5.5

    elif month == 3: #MARCH
        H_0 = 33.86
        A = 0.5486
        B = -0.1996
        C = 0.0280
        D = 0.0216
        E = -0.0176
        epsilon = 9.4
        F = -0.04
        sigma_f = 0.10
        IAV = 6.2

    elif month == 4: #APRIL
        H_0 = 36.80
        A = 0.4046
        B = -0.4005
        C = 0.0256
        D = 0.0024
        E = -0.0641
        epsilon = 9.4
        F = -0.09
        sigma_f = 0.09
        IAV = 6.1

    elif month == 5: #MAY
        H_0 = 36.93
        A = 0.0214
        B = -1.1795
        C = -0.1076
        D = -0.0244
        E = -0.0142
        epsilon = 10.6
        F = -0.21
        sigma_f = 0.09
        IAV = 6.3

    elif month == 6: #JUNE
        H_0 = 36.59
        A = 0.7021
        B = -1.4819
        C = -0.1195
        D = -0.0009
        E = -0.0603
        epsilon = 14.1
        F = -0.16
        sigma_f = 0.12
        IAV = 8.1
 
    elif month == 7: # JULY
        H_0 = 11.02
        A = 0.3008
        B = -1.2591
        C = -0.0811
        D = -0.0043
        E = -0.0959
        epsilon = 9.5
        F = 0.02
        sigma_f = 0.10
        IAV = 6.7

    elif month == 8: # AUGUST
        H_0 = 4.64
        A = 0.3100
        B = -0.6350
        C = -0.0655
        D = 0.0059
        E = -0.0005
        epsilon = 4.6
        F = -0.01
        sigma_f = 0.05
        IAV = 3.3

    elif month == 9: # SEPTEMBER
        H_0 = 15.81
        A = 0.2119
        B = -1.0292
        C = -0.0868
        D = -0.0177
        E = -0.0723
        epsilon = 7.8
        F = -0.03
        sigma_f = 0.06
        IAV = 3.8

    elif month == 10: # OCTOBER
        H_0 = 22.66
        A = 0.3594
        B = -1.3483
        C = -0.1063
        D = 0.0051
        E = -0.0577
        epsilon = 8.0
        F = -0.08
        sigma_f = 0.06
        IAV = 4.0

    elif month == 11: # NOVEMBER
        H_0 = 25.57
        A = 0.1496
        B = -1.4643
        C = -0.1409
        D = -0.0079
        E = -0.0258
        epsilon = 7.9
        F = -0.05
        sigma_f = 0.07
        IAV = 4.3

    elif month == 12: #DECEMBER
        H_0 = 26.67
        A = -0.1876
        B = -1.4229
        C = -0.1413
        D = -0.0316
        E = -0.0029
        epsilon = 8.2
        F = -0.06
        sigma_f = 0.07
        IAV = 4.8

    else:
        raise ValueError('Not a valid integer for a month!')

    snow_depth = H_0 + (A*X) + (B*Y) + (C*X*Y) + (D*X*X) + (E*Y*Y)

    return snow_depth, epsilon

def SWE(lat, lon, month):

    from numpy import sin, cos, pi

    Y = (90-lat) * sin(lon*(pi/180))
    X = (90-lat) * cos(lon*(pi/180))


    #Set the coefficients

    if month == 1: #JANUARY
        H_0 = 8.57
        A = -0.0270
        B = -0.3400
        C = -0.0319
        D = -0.0056
        E = -0.0005
        epsilon = 2.5
        F = -0.05
        sigma_f = 0.024
        IAV = 1.6

    elif month == 2: #FEBRUARY
        H_0 = 9.45
        A = 0.0058
        B = -0.1309
        C = 0.0017
        D = -0.0021
        E = -0.0072
        epsilon = 2.6
        F = -0.07
        sigma_f = 0.028
        IAV = 1.8

    elif month == 3: #MARCH
        H_0 = 10.74
        A = 0.1618
        B = 0.0276
        C = 0.0213
        D = 0.0076
        E = -0.0125
        epsilon = 3.1
        F = 0.07
        sigma_f = 0.032
        IAV = 2.1

    elif month == 4: #APRIL
        H_0 = 11.67
        A = 0.0841
        B = -0.1328
        C = 0.0081
        D = -0.0003
        E = -0.0301
        epsilon = 3.2
        F = -0.013
        sigma_f = 0.032
        IAV = 2.1

    elif month == 5: #MAY
        H_0 = 11.80
        A = -0.0043
        B = -0.4284
        C = -0.0380
        D = -0.0071
        E = -0.0063
        epsilon = 3.5
        F = -0.047
        sigma_f = 0.033
        IAV = 2.2

    elif month == 6: #JUNE
        H_0 = 12.48
        A = 0.2084
        B = -0.5739
        C = -0.0468
        D = -0.0023
        E = -0.0253
        epsilon = 4.9
        F = -0.030
        sigma_f = 0.044
        IAV = 2.9
 
    elif month == 7: # JULY
        H_0 = 4.01
        A = 0.0970
        B = -0.4930
        C = -0.0333
        D = -0.0026
        E = -0.0343
        epsilon = 3.5
        F = 0.008
        sigma_f = 0.037
        IAV = 2.4

    elif month == 8: # AUGUST
        H_0 = 1.08
        A = 0.0712
        B = -0.1450
        C = -0.0155
        D = 0.0014
        E = -0.0000
        epsilon = 1.1
        F = -0.001
        sigma_f = 0.012
        IAV = 0.8

    elif month == 9: # SEPTEMBER
        H_0 = 3.84
        A = 0.0393
        B = -0.2107
        C = -0.0182
        D = -0.0053
        E = -0.0190
        epsilon = 2.0
        F = -0.003
        sigma_f = 0.016
        IAV = 1.0

    elif month == 10: # OCTOBER
        H_0 = 6.24
        A = 0.1158
        B = -0.2803
        C = -0.0215
        D = 0.0015
        E = -0.0176
        epsilon = 2.3
        F = -0.005
        sigma_f = 0.021
        IAV = 1.4

    elif month == 11: # NOVEMBER
        H_0 = 7.54
        A = 0.0567
        B = -0.3201
        C = -0.0284
        D = -0.0032
        E = -0.0129
        epsilon = 2.4
        F = -0.000
        sigma_f = 0.023
        IAV = 1.5

    elif month == 12: #DECEMBER
        H_0 = 8.00
        A = -0.0540
        B = -0.3650
        C = -0.0362
        D = -0.0112
        E = -0.0035
        epsilon = 2.5
        F = -0.003
        sigma_f = 0.024
        IAV = 1.5

    else:
        raise ValueError('Not a valid integer for a month!')

    snow_depth = H_0 + (A*X) + (B*Y) + (C*X*Y) + (D*X*X) + (E*Y*Y)

    return snow_depth, epsilon
