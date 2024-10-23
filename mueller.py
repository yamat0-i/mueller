import numpy as np
import matplotlib.pyplot as plt


def cos(x):
    if x == np.pi/2 or x == (3/2)*np.pi:
        y = 0
    else:
        y = np.cos(x)
    return y

def sin(x):
    if x == 0 or x == np.pi:
        y = 0
    else:
        y = np.sin(x)
    return y

def get_stokes_vector(I, phi, chi):
    """Calc Stokes vector: S = (s0, s1, s2, s3)

    var:
        I: intensity of light. If normalized, I=1.
        phi: longitude angle
        chi: latitude angle
    """
    S = np.array([
        [I],
        [I * cos(2*chi) * cos(2*phi)],
        [I * cos(2*chi) * sin(2*phi)],
        [I * sin(2*chi)],
    ])

    return S

def get_mueller_matrix(theta, ret):
    """Calc Mueller matrix.

    var:
        theta: Fast axes angle of waveplate
        ret: Retardance
    """
    M = np.matrix([
        [1, 0, 0, 0],
        [0, cos(2*theta)**2 + cos(ret)*sin(2*theta)**2, (1 - cos(ret))*sin(2*theta)*cos(2*theta), sin(ret)*sin(2*theta)],
        [0, (1 - cos(ret))*sin(2*theta)*cos(2*theta), sin(2*theta)**2 + cos(ret)*cos(2*theta)**2, -sin(ret)*cos(2*theta)],
        [0, -sin(ret)*sin(2*theta), sin(ret)*cos(2*theta), cos(ret)]])
    
    return M

if __name__ == '__main__':
    # S0: Input pol
    # H state: (1, 0, 0)
    # V state: (1, np.pi/2, 0)
    # D state: (1, np.pi/4, 0)
    # A state: (1, -np.pi/4, 0)
    # R state: (1, 0, np.pi/4)
    # L state: (1, 0, -np.pi/4)
    S0  = get_stokes_vector(1, 0, 0)
    print('S0', S0)

    # Retardance of waveplates
    retardance1 = np.pi/2
    retardance2 = np.pi/2
    retardance3 = 0

    # Range of waveplate's angle
    theta1 = np.linspace(0, 2*np.pi, 100)
    theta2 = np.linspace(0, 2*np.pi, 100)
    theta3 = np.linspace(0, np.pi, 90)

    # Mueller matrix
    # M_QWP = get_mueller_matrix(theta, np.pi/2)
    # M_HWP = get_mueller_matrix(theta, np.pi)

    # Calc output pol
    # M2 <- M1 <- Hstate
    # S2 = M2(S1) = M2(M1(S0)): Stokes vector, 'list' of matrix
    # s1, s2, s3: Stokes parameters, 'list'
    print('Calculating...')
    S2 = []
    s1 = []
    s2 = []
    s3 = []
    for i in range(len(theta2)):
        for j in range(len(theta1)):
            M1 = get_mueller_matrix(theta1[j], retardance1)
            S1_value = np.dot(M1, S0)    
            M2 = get_mueller_matrix(theta2[i], retardance2)
            S2_value = np.dot(M2, S1_value)
            S2.append(S2_value)
            s1.append(S2_value[1,0])
            s2.append(S2_value[2,0])
            s3.append(S2_value[3,0])


    # S3 = []
    # for i in range(len(theta2)):
    #     for j in range(len(theta1)):
    #         for k in range(len(theta3)):
    #             M1 = get_mueller_matrix(theta3[k], np.pi/2)
    #             S1_value = np.dot(M1, S0)
    #             M2 = get_mueller_matrix(theta1[j], np.pi/2)
    #             S2_value = np.dot(M2, S1_value)    
    #             M3 = get_mueller_matrix(theta3[i], np.pi)
    #             S3_value = np.dot(M3, S2_value)
    #             S3.append(S3_value)
    #             s1.append(S3_value[1,0])
    #             s2.append(S3_value[2,0])
    #             s3.append(S3_value[3,0])
            
    N = len(s1) # s1, s2, s3 are arrays of Stokes parameters obtained in LTM
    psi = np.empty(N)
    kai = np.empty(N)
    for i in range(N):
        s12 = s1[i]**2 + s2[i]**2
        psi[i] = (1/2) * np.arctan2(s2[i], s1[i]) * (180/np.pi)
        kai[i] = (1/2) * np.arctan2(s3[i], np.sqrt(s12)) * (180/np.pi)
        np.append(psi, psi[i])
        np.append(kai, kai[i])


    print('Creating a plot...')
    # Plot
    # Poincare sphere
    r = 1
    theta_1_0 = np.linspace(0, 2*np.pi, 400)
    theta_2_0 = np.linspace(0, 2*np.pi, 400)
    theta_1, theta_2 = np.meshgrid(theta_1_0, theta_2_0)
    x = np.cos(theta_2)*np.sin(theta_1) * r
    y = np.sin(theta_2)*np.sin(theta_1) * r
    z = np.cos(theta_1) * r

    fig = plt.figure(layout='tight')
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(x,y,z, alpha=0.05)
    plt.xlim([1,-1])
    plt.ylim([1,-1])
    ax.set_zlim([-1,1])
    ax.set_aspect('equal')

    # Grid
    x_cs2 = np.cos(0)*np.sin(theta_1_0) * r
    y_cs2 = np.sin(0)*np.sin(theta_1_0) * r
    z_cs2 = np.cos(theta_1_0) * r

    ax.plot(x_cs2, y_cs2, z_cs2, color="black", linewidth=0.5)

    x_cs1 = np.cos(np.pi/2)*np.sin(theta_1_0) * r
    y_cs1 = np.sin(np.pi/2)*np.sin(theta_1_0) * r
    z_cs1 = np.cos(theta_1_0) * r

    ax.plot(x_cs1, y_cs1, z_cs1, color="black", linewidth=0.5)

    x_cs3 = np.cos(theta_2_0)*np.sin(np.pi/2) * r
    y_cs3 = np.sin(theta_2_0)*np.sin(np.pi/2) * r
    z_cs3 = np.cos(np.pi/2) * r

    ax.plot(x_cs3, y_cs3, z_cs3, color="black", linewidth=0.5)

    # data
    ax.scatter(s1, s2, s3, marker='.', color='blue', label='10x', alpha=0.5)

    ax.set_xlabel('s1', fontsize=18)
    ax.set_ylabel('s2', fontsize=18)
    ax.set_zlabel('s3', fontsize=18)

    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='z', labelsize=8)

    # ax.set_title(r'H$\rightarrow$QWP$\rightarrow$QWP')

    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.scatter(psi, kai, marker='o', alpha=0.1)
    ax1.set_xlabel(r'$\psi$ [deg]')
    ax1.set_ylabel(r'$\chi$ [deg]')

    print('Finished.')
    plt.show(block=False)
    input()
