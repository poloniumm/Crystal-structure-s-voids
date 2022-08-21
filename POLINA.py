import numpy
import math
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if len(sys.argv) == 15:
    cif = open(sys.argv[12], 'r')
    f = [line.rstrip() for line in cif] #Здесь начинаем считывание данных из cif-файла
    formula = f[0]
    #Параметры ячейки, углы(находятся в 4 строке)
    parametrs = f[3].split()
    a, b, c, alpha, beta, gamma = float(parametrs[0]), float(parametrs[1]), float(parametrs[2]), float(parametrs[3]), float(parametrs[4]), float(parametrs[5])
    #ФУНКЦИЯ РАССТОЯНИЙ
    def distance (x, y, z, a, b, c, m1, l1, n1, alpha, beta, gamma):
        dist = math.sqrt(((x - l1) * a) ** 2 + ((y - m1) * b) ** 2 + ((z - n1) * c) ** 2 + 2 * ((x - l1) * a) * (
                    (y - m1) * b) * math.cos(math.radians(gamma)) + 2 * ((x - l1) * a) * ((z - n1) * c) * math.cos(
            math.radians(beta)) + 2 * ((y - m1) * b) * ((z - n1) * c) * math.cos(math.radians(alpha)))
        return dist

    def sefe (xse1, yse1, zse1, xse2, yse2, zse2, alpha, beta, gamma):
        dsefe = ((xse1 - xse2) ** 2 + (yse1 - yse2) ** 2 + (zse1 - zse2) ** 2 + (zse1 - zse2) * (
                       xse1 - xse2) * math.cos(math.radians(beta)) + (zse1 - zse2) * (yse1 - yse2) * math.cos(
              math.radians(alpha)) + (xse1 - xse2) * (yse1 - yse2) * math.cos(math.radians(gamma))) ** 0.5
        return dsefe
    #Создаём список атомов и их координат по отдельным строкам
    odinatom = []
    atoms = []
    for i in range(6, len(f)):
        atom = []
        atom.append(f[i].split())
        atom = atom[0]
        atoms.append(atom)
    #Определяем начало и конец отсчёта (относительно)
    Xstart = float(sys.argv[1])
    Ystart = float(sys.argv[2])
    Zstart = float(sys.argv[3])
    Xend = float(sys.argv[4])
    Yend = float(sys.argv[5])
    Zend = float(sys.argv[6])
    step = float(sys.argv[7])
    start1 = float(sys.argv[8])
    end1 = str(sys.argv[9])
    start2 = str(sys.argv[10])
    end2 = str(sys.argv[11])
    reswith = open(sys.argv[13], 'w')
    reswithout = open(sys.argv[14], 'w')
    #Делаем шапку с ограничениями в res - файлах
    reswithout.write(formula + '- test without radii' + '\n' + '\n')
    reswithout.write("Стартовая точка исследования (в долях):" + str(Xstart) + ", " + str(Ystart) + ", " + str(Zstart) + '\n')
    reswithout.write("Конечная точка исследования (в долях):" + str(Xend) + ", " + str(Yend) + ", " + str(Zend) + '\n' + '\n' + '\n')
    reswith.write(formula + " - test with radius" + '\n' + "\n")
    reswith.write("Стартовая точка исследования (в долях):" + str(Xstart) + ", " + str(Ystart) + ", " + str(Zstart) + '\n')
    reswith.write("Конечная точка исследования (в долях):" + str(Xend) + ", " + str(Yend) + ", " + str(Zend) + '\n' + '\n' + '\n')
    if str(sys.argv[9]) == '~':
        reswithout.write('Condition: d>' + str(start1) + '\n')
    else:
        reswithout.write('Condition: ' + str(start1) + '<d<' + str(end1) + '\n')
    if str(sys.argv[11]) == '~':
        reswith.write('Condition: d>' + str(start2) + '\n')
    else:
        reswith.write('Condition: ' + str(start2) + '<d<' + str(end2) + '\n')
#Расчет расстояний
    fi = []
    se = []
    l1 = Xstart
    while l1 <= Xend:
        m1 = Ystart
        while m1 <= Yend:
            n1 = Zstart
            while n1 <= Zend:
                schetchik1 = 0
                schetchik2 = 0
                sd = []
                for i in atoms:
                    lable = i[0]

                    x, y, z, radius = float(i[1]), float(i[2]), float(i[3]), float(i[4])
                    if l1 - x >= 0.5: x += 1
                    elif l1 - x < -0.5: x -= 1
                    if m1 - y > 0.5: y += 1
                    elif m1 - y < -0.5: y -= 1
                    if n1 - z > 0.5: z += 1
                    elif n1 - z < -0.5: z -= 1
                    #вводим формулу расстояний
                    #dist = math.sqrt(((x - l1)*a) ** 2 + ((y - m1)*b) ** 2 + ((z - n1)*c) ** 2 + 2 * ((x - l1)*a) * ((y - m1)*b) * math.cos(math.radians(gamma)) + 2 * ((x - l1)*a) * ((z - n1)*c) * math.cos(math.radians(beta)) + 2 * ((y - m1)*b) * ((z - n1)*c) * math.cos(math.radians(alpha)))


                    sd.append([lable, distance(x, y, z, a, b, c, m1, l1, n1, alpha, beta, gamma)])
                    # 1случай - [x,беск)
                    if str(sys.argv[9]) == "~":
                        end1 = str(sys.argv[9])
                        if start1 < distance(x, y, z, a, b, c, m1, l1, n1, alpha, beta, gamma):
                           schetchik1 += 1
                        if schetchik1 == len(atoms):
                           l1 = float('{:.4f}'.format(l1))
                           m1 = float('{:.4f}'.format(m1))
                           n1 = float('{:.4f}'.format(n1))

                           reswithout.write("(" + str(l1) + ", " + str(m1) + ", " + str(n1) + ')' + '\n')
                           for k in range(len(atoms)):
                             reswithout.write('dist from point to atom '+str(atoms[k][0])+' '+ ' [' + str(atoms[k][1]) +',' +str(atoms[k][2])+','+str(atoms[k][3])+ '] ' + " = " + str(sd[k][1]) + '\n')
                           fi.append((float(l1)*a, float(m1)*b, float(n1)*c))
                           schetchik1 = 0
                    else:
                       # 2 случай - (0,x]
                        end1 = float(sys.argv[9])
                        if start1 < dist < end1:
                           schetchik1 += 1
                        if schetchik1 == len(atoms):
                           l1 = float('{:.4f}'.format(l1))
                           m1 = float('{:.4f}'.format(m1))
                           n1 = float('{:.4f}'.format(n1))
                           reswithout.write("(" + str(l1) + ", " + str(m1) + ", " + str(n1) + ')' + '\n')
                           for k in range(len(atoms)):
                                   reswithout.write('dist from point to atom ' +str(atoms[k][0])+' '+ ' [' + str(atoms[k][1]) +',' +str(atoms[k][2])+','+str(atoms[k][3])+ '] ' + " = " + str(sd[k][1]) + '\n')
                           schetchik1 = 0
                           fi.append((float(l1)*a, float(m1)*b, float(n1)*c))
                    end2 = float(sys.argv[10])
                # 1 случай б.у.
                    if str(sys.argv[11]) == '~':
                        end2 = str(sys.argv[11])
                        if float(start2) < (float(distance(x, y, z, a, b, c, m1, l1, n1, alpha, beta, gamma)) - float(radius)):
                            schetchik2 += 1
                        if schetchik2 == len(atoms):
                            l1 = float('{:.4f}'.format(l1))
                            m1 = float('{:.4f}'.format(m1))
                            n1 = float('{:.4f}'.format(n1))
                            reswith.write('(' + str(l1) + ", " + str(m1) + ", " + str(n1) + ')' + '\n')
                            for k in range(len(atoms)):
                                   reswith.write('dist from point to atom ' + str(atoms[k][0])+' '+ ' [' + str(atoms[k][1]) +',' +str(atoms[k][2])+','+str(atoms[k][3])+ '] ' + " = " + str(sd[k][1]) + '\n')
                            se.append((float(l1)*a, float(m1)*b, float(n1)*c))
                            schetchik2 = 0
                    else:
                        end2 = float(sys.argv[11])
                        if float(start2) < float(distance(x, y, z, a, b, c, m1, l1, n1, alpha, beta, gamma)) - float(radius) < end2:
                            schetchik2 += 1
                        if schetchik2 == len(atoms):
                            l1 = float('{:.4f}'.format(l1))
                            m1 = float('{:.4f}'.format(m1))
                            n1 = float('{:.4f}'.format(n1))
                            reswith.write('(' + str(l1) + ", " + str(m1) + ", " + str(n1) + ')' + '\n')
                            for k in range(len(atoms)):
                                   reswith.write('dist from point to atom '+str(atoms[k][0])+' ' + ' [' + str(atoms[k][1]) +',' +str(atoms[k][2])+','+str(atoms[k][3])+ '] ' + " = " + str(sd[k][1]) + '\n')
                            se.append((float(l1)*a, float(m1)*b, float(n1)*c))
                            schetchik2 = 0
                n1 += step
            m1 += step
        l1 += step
        m1 = float(sys.argv[2])
        n1 = float(sys.argv[3])
    cif.close()
    reswithout.close()
    reswith.close()
    print('Session ended (' + sys.argv[13] + ', ' + sys.argv[14] + ')')

    atoms_print = []
    for i in range (len(atoms)):
        name_ = atoms[i][0]
        x_ = round(float(atoms[i][1]) * float(a), 3)
        y_ = round(float(atoms[i][2]) * float(b), 3)
        z_ = round(float(atoms[i][3]) * float(c), 3)
        ar_ = atoms[i][4]
        at_= [name_, x_, y_, z_, ar_]
        atoms_print.append(at_)


    # определение соседей для подошедших точек
    pairs_se = []
    pairs_fi = []


    # for i in range(len(se)):
    #     xse1 = se[i][0]
    #     yse1 = se[i][1]
    #     zse1 = se[i][2]
    #     for k in range(len(se)):
    #         xse2 = se[k][0]
    #         yse2 = se[k][1]
    #         zse2 = se[k][2]
    #
    #         distse = ((xse1-xse2)**2 + (yse1-yse2)**2 + (zse1-zse2)**2 + (zse1-zse2)*(xse1-xse2)*math.cos(math.radians(beta)) + (zse1-zse2)*(yse1-yse2)*math.cos(math.radians(alpha)) + (xse1-xse2)*(yse1-yse2)*math.cos(math.radians(gamma)))**0.5
    #
    #         if distse < float(step)*((a+b+c)/3)*1.6 and distse > float(step)*((a+b+c)/3)*0.2:
    #             pair = [se[i], se[k]]
    #             pairs_se.append(pair)


    # for i in range(len(fi)):
    #     xfi1 = fi[i][0]
    #     yfi1 = fi[i][1]
    #     zfi1 = fi[i][2]
    #     for k in range (len(fi)):
    #         xfi2 = fi[k][0]
    #         yfi2 = fi[k][1]
    #         zfi2 = fi[k][2]
    #
    #         distfi = ((xfi1-xfi2)**2 + (yfi1-yfi2)**2 + (zfi1-zfi2)**2 + (zfi1-zfi2)*(xfi1-xfi2)*math.cos(math.radians(beta)) + (zfi1-zfi2)*(yfi1-yfi2)*math.cos(math.radians(alpha)) + (xfi1-xfi2)*(yfi1-yfi2)*math.cos(math.radians(gamma)))**0.5
    #         if distfi < float(step)*((a+b+c)/3)*1.6 and distfi > float(step)*((a+b+c)/3)*0.2:
    #             pair = [fi[i], fi[k]]
    #             pairs_fi.append(pair)


    # визуализация
    try:
        fig = plt.figure()
        ax = fig.add_subplot(121, projection='3d')
        [atn, xa, ya, za, ar] = list(map(list, zip(*atoms_print)))
        for i, item in enumerate(xa):
            xa[i] = float(item)
        for i, item in enumerate(ya):
            ya[i] = float(item)
        for i, item in enumerate(za):
            za[i] = float(item)
        ax.scatter(xa, ya, za, color='purple')
        for i in atoms_print:
            ax.text(float(i[1]), float(i[2]), float(i[3]), i[0])
        [x, y, z] = list(map(list, zip(*fi)))
        l = ax.scatter(x, y, z)

        Xfi = []
        Yfi = []
        Zfi = []
        for i in range(len(fi)):
            xfi1 = fi[i][0]
            yfi1 = fi[i][1]
            zfi1 = fi[i][2]
            for k in range(len(fi)):
                xfi2 = fi[k][0]
                yfi2 = fi[k][1]
                zfi2 = fi[k][2]

                #distfi = ((xfi1 - xfi2) ** 2 + (yfi1 - yfi2) ** 2 + (zfi1 - zfi2) ** 2 + (zfi1 - zfi2) * (
                          #  xfi1 - xfi2) * math.cos(math.radians(beta)) + (zfi1 - zfi2) * (yfi1 - yfi2) * math.cos(
                    #math.radians(alpha)) + (xfi1 - xfi2) * (yfi1 - yfi2) * math.cos(math.radians(gamma))) ** 0.5


                if float(sefe(xfi1, yfi1, zfi1, xfi2, yfi2, zfi2, alpha, beta, gamma)) < (float(step) * ((a**2 + b**2 + c**2+ a*b*math.cos(math.radians(gamma) + a*c*math.cos(math.radians(beta)
                                  +b*c*math.cos(math.radians(gamma))))) ** 0.5) * 1.1) and sefe(xfi1, yfi1, zfi1, xfi2, yfi2, zfi2, alpha, beta, gamma) > ((float(step) * (a**2 + b**2 + c**2
                                + a*b*math.cos(math.radians(gamma)) + a*c*math.cos(math.radians(beta))
                                  +b*c*math.cos(math.radians(gamma))) ** 0.5) * 0.01):
                    x, y, z = (xfi1, xfi2), (yfi1, yfi2), (zfi1, zfi2)
                    ax.plot(x, y, z, color='r')


    except:
        print("No one for the first calculation")
    try:
        ax = fig.add_subplot(122, projection='3d')
        [atn, xa, ya, za, ar] = list(map(list, zip(*atoms_print)))
        for i, item in enumerate(xa):
            xa[i] = float(item)
        for i, item in enumerate(ya):
            ya[i] = float(item)
        for i, item in enumerate(za):
            za[i] = float(item)
        ax.scatter(xa, ya, za, color='purple')
        for i in atoms_print:
            ax.text(float(i[1]), float(i[2]), float(i[3]), i[0])
        [x, y, z] = list(map(list, zip(*se)))
        l = ax.scatter(x, y, z)


        Xse = []
        Yse = []
        Zse = []
        for i in range(len(se)):
            xse1 = se[i][0]
            yse1 = se[i][1]
            zse1 = se[i][2]
            for k in range(len(se)):
                xse2 = se[k][0]
                yse2 = se[k][1]
                zse2 = se[k][2]
                r1=1
                r2=1
                r3=1
                #distse = ((xse1 - xse2) ** 2 + (yse1 - yse2) ** 2 + (zse1 - zse2) ** 2 + (zse1 - zse2) * (
                 #           xse1 - xse2) * math.cos(math.radians(beta)) + (zse1 - zse2) * (yse1 - yse2) * math.cos(
                  #  math.radians(alpha)) + (xse1 - xse2) * (yse1 - yse2) * math.cos(math.radians(gamma))) ** 0.5
                #
                if float(sefe(xse1, yse1, zse1, xse2, yse2, zse2, alpha, beta, gamma)) < (float(step) * ((a ** 2 + b ** 2 + c ** 2 + a * b * math.cos(math.radians(gamma)
                 + a * c * math.cos(math.radians(beta)+ b * c * math.cos(math.radians(gamma))))) ** 0.5) * 1.1) and sefe(xse1, yse1, zse1, xse2, yse2, zse2, alpha, beta, gamma)  > ((float(step) * (a ** 2 + b ** 2 + c ** 2+ a * b * math.cos(math.radians(gamma)) + a * c * math.cos(math.radians(beta))
                             + b * c * math.cos(math.radians(gamma))) ** 0.5) * 0.01):
                    x, y, z = (xse1, xse2), (yse1, yse2), (zse1, zse2)
                    ax.plot(x, y, z, color='r')

    except:
        print("No one for the second calculation")
    plt.show()
else:
    print("Error")
