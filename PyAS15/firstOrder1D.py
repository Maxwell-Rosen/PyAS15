''' Maxwell Rosen November 27 2023
This file contains the implementation of the IAS15 algorithm for the first order ODE of the form
dy/dt = f(y).
The algorithm is implemented in the function pushIAS15(derivative, startingPoint, dt, lastB=None)
The goal is to eventually extend this to a PDE, so we will first solve a simpler ODE system.
One test is to integrate some sample ODE's with these functions.
Examples: Exponenetial growth, simple harmonic oscillator, logistic growth
Ultimately, aim for the equations
dy/ds = d psi(x,y)/ dx = f(x,y)
dx/ds = -d psi(x,y)/ dy= g(x,y)
'''
import numpy as np

global h, rr, c, d, w

h = np.array([0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648,\
                0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626])

rr = np.array([0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278, 0.3526247171131696373739078, \
                0.2963621565762474909082556, 0.1723840253762772723863278, 0.5471536263305553830014486, 0.4908910657936332365357964, \
                0.3669129345936630180138686, 0.1945289092173857456275408, 0.7342101772154105315232106, 0.6779476166784883850575584, \
                0.5539694854785181665356307, 0.3815854601022408941493028, 0.1870565508848551485217621, 0.8853209468390957680903598, \
                0.8290583863021736216247076, 0.7050802551022034031027798, 0.5326962297259261307164520, 0.3381673205085403850889112, \
                0.1511107696236852365671492, 0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035945, \
                0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769703679639, 0.0921996667221917338008147])

c = np.array([-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, -0.0035758977292516175949345, \
                0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006, -0.0547553868890686864408084, \
                0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915424460, 0.0421585277212687077072973, \
                -0.3600995965020568122897665, 1.2501507118406910258505441, -1.8704917729329500633517991, 0.0012717903090268677492943,\
                -0.0387603579159067703699046, 0.3609622434528459832253398, -1.4668842084004269643701553, 2.9061362593084293014237913, \
                -2.7558127197720458314421588])

d = np.array([0.0562625605369221464656522, 0.0031654757181708292499905, 0.2365032522738145114532321, 0.0001780977692217433881125,\
                0.0457929855060279188954539, 0.5891279693869841488271399, 0.0000100202365223291272096, 0.0084318571535257015445000, \
                0.2535340690545692665214616, 1.1362815957175395318285885, 0.0000005637641639318207610, 0.0015297840025004658189490, \
                0.0978342365324440053653648, 0.8752546646840910912297246, 1.8704917729329500633517991, 0.0000000317188154017613665,\
                0.0002762930909826476593130, 0.0360285539837364596003871, 0.5767330002770787313544596, 2.2485887607691597933926895, \
                2.7558127197720458314421588])

w = np.array([0.03125, 0.185358154802979278540728972807180754479812609, 0.304130620646785128975743291458180383736715043, \
                0.376517545389118556572129261157225608762708603, 0.391572167452493593082499533303669362149363727, 0.347014795634501068709955597003528601733139176,\
                0.249647901329864963257869294715235590174262844, 0.114508814744257199342353731044292225247093225])

def push(derivative: callable, startingPoint, dt, B = np.zeros(7)):
    ''' Input:
    derivative - the derivative of the function to integrate. Right hand side of the ODE
    startingPoint - the starting point of the integration. x_final
    Output:
    newPosition - the new position after the integration. x_final'''

    for i in range(12):
        xh = calculateNodePositions(B, startingPoint, derivative, dt, hp = h)
        F = calculateDerivatives(derivative, xh)
        G = calculateGFromF(F)
        B = calculateB(G)
        #print("B", B, "G", G, "F", F, "xh", xh)
       #print("B", B)
    newPosition = calculateNewPosition(B, startingPoint, derivative, dt)
    return newPosition, B

def calculateNodePositions(B, startingPoint, derivative: callable, dt, hp = h):
    F1 = derivative(startingPoint)
    xh = np.zeros(len(hp))
    for n in range(len(hp)):
        xh[n] = startingPoint + hp[n]*dt*(F1   + hp[n]*1/2*\
                                        ( B[0] + hp[n]*2/3*\
                                        ( B[1] + hp[n]*3/4*\
                                        ( B[2] + hp[n]*4/5*\
                                        ( B[3] + hp[n]*5/6*\
                                        ( B[4] + hp[n]*6/7*\
                                        ( B[5] + hp[n]*7/8*\
                                          B[6])))))))
    return xh

def calculateNewPosition(B, startingPoint, derivative: callable, dt):
    newPosition = calculateNodePositions(B, startingPoint, derivative, dt, hp = np.array([1]))
    return newPosition

def calculateDerivatives(derivative, x):
    F = np.zeros(len(x))
    for i in range(len(x)):
        F[i] = derivative(x[i])
    return F


def calculateB(G):
    B = np.zeros(len(G))
    B[0] = G[0] + c[0]*G[1] + c[1]*G[2] + c[3]*G[3] + c[6]*G[4] + c[10]*G[5] + c[15]*G[6]
    B[1] =             G[1] + c[2]*G[2] + c[4]*G[3] + c[7]*G[4] + c[11]*G[5] + c[16]*G[6]
    B[2] =                         G[2] + c[5]*G[3] + c[8]*G[4] + c[12]*G[5] + c[17]*G[6]
    B[3] =                                     G[3] + c[9]*G[4] + c[13]*G[5] + c[18]*G[6]
    B[4] =                                                 G[4] + c[14]*G[5] + c[19]*G[6]
    B[5] =                                                              G[5] + c[20]*G[6]
    B[6] =                                                                           G[6]
    return B

def calculateGFromF(F):
    G = np.zeros(len(F)-1)
    G[0] =       (F[1] - F[0])/rr[0]
    G[1] =      ((F[2] - F[0])/rr[1]  - G[0])/rr[2]
    G[2] =     (((F[3] - F[0])/rr[3]  - G[0])/rr[4]  - G[1])/rr[5]
    G[3] =    ((((F[4] - F[0])/rr[6]  - G[0])/rr[7]  - G[1])/rr[8]  - G[2])/rr[9]
    G[4] =   (((((F[5] - F[0])/rr[10] - G[0])/rr[11] - G[1])/rr[12] - G[2])/rr[13] - G[3])/rr[14]
    G[5] =  ((((((F[6] - F[0])/rr[15] - G[0])/rr[16] - G[1])/rr[17] - G[2])/rr[18] - G[3])/rr[19] - G[4])/rr[20]
    G[6] = (((((((F[7] - F[0])/rr[21] - G[0])/rr[22] - G[1])/rr[23] - G[2])/rr[24] - G[3])/rr[25] - G[4])/rr[26] - G[5])/rr[27]
    return G