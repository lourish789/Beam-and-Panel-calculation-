
import math
import numpy as np
class raft_design:
    def __init__(self):
        pass
    def panel(self):
            y = int(input("How many panels are involved in the structure "))
            e = round(float(input("What is the eccentricity ")), 2)
            d = round(float(input("What is the base diameter of the column size ")), 2)
            b = round(float(input("What is the base depth of the column size ")), 2)
            w = round(float(input("What is the total load ")), 2)
            a = b * d  #area
            P1 = round((w / a) + ((6 * w * e) / (b * (d ** 2))), 2)
            P2 = round((w / a) - ((6 * w * e) / (b * (d ** 2))), 2)
            print("The Pressure on the Column size is ", P1, " or ", P2)
            print("okay")
            # at Ultimate limit safety
            P1 = round((P1 * 1.47), 2)
            P2 = round(float(P2 * 1.47), 2)
            # Pressure at gridline(Px)
            aa = round(float(input("What is the distance between column a and b ")), 2)
            bb = round(float(input("What is the distance between column b and c ")), 2)
            cc = round(float(input("What is the distance between column c and d ")), 2)
            abc = round(aa + bb + cc, 2)
            # Pressure value at grid line using linear principle
            Pa = round(((aa / abc) * (P2 - P1) + P1), 2)
            Pb = round(float((bb / abc) * (P2 - P1) + Pa), 2)
            Pc = round(float((cc / abc) * (P2 - P1) + Pb), 2)
            plist = [P1, Pa, Pb, Pc]
            pdict = {
                "Grid A": P1,
                "Grid B": Pa,
                "Grid C": Pb,
                "Grid D": Pc
            }
            for p in pdict.keys():
                print("The Pressure at ", p, " is ", + pdict[p])
                # print(plist)
            # concrete own pressure
            cop = round(0.2 * 24 * 1.4, 2)
            print("Assuming the slab is 200mm thick, the concrete own pressure will be: ", cop)
            # panel()
            y = int(input("How many panels are involved in the structure  "))
            # x = range(1,y+1)
            for i in range(y):
                Answers = {}
                if i <= y:
                    f = int(input(
                        "If this panel is designed as a two-way slab with UDL, enter 1, else enter any other number: "))
                    if f == 1:
                        print('okay, ')
                        ly = round(float(input("For the this panel, What is the value of ly; ")), 2)
                        lx = round(float(input("what is the value of lx ")), 2)
                        if ly / lx < 2:
                            s1 = round(
                                float(input("Input the coefficient of the continuous edge of the shorter span: ")), 2)
                            s2 = round(float(input("Input the coefficient of the mid span of the shorter span: ")), 2)
                            s3 = round(
                                float(input("Input the coefficient of the continuous edge of the longer span: ")), 2)
                            s4 = round(float(input("Input the coefficient of the mid span of the longer span:  ")), 2)
                            # shorter span
                            # midspan
                            d_l = round(pdict[p] - cop, 2)  # design load(pressure)
                            Fcu = 25
                            b = 1000
                            m1 = round(s2 * d_l * (lx ** 2), 2)
                            k = round(m1 / ((Fcu * b) * d ** 2), 2)
                            La = round(0.5 + math.sqrt(0.25 - (k / 0.9)), 2)
                            Z = round(La * d, 2)
                            Fy = 410
                            Ans1 = m1 / (0.95 * Fy * Z)
                            Ans1 = ("The area of the mid span for the first panel under the shorter span is: ", + Ans1)
                            print(Ans1)
                            # continuos edge
                            m2 = round(s1 * d_l * (lx ** 2), 2)
                            k = round(m2 / ((Fcu * b) * d ** 2), 2)
                            La = round(0.5 + math.sqrt(0.25 - (k / 0.9)), 2)
                            Z = La * d
                            Ans2 = m2 / (0.95 * Fy * Z)
                            Ans2 = (
                            "Th area of the continuous edge for the first panel under the shorter span is: ", + Ans2)
                            print(Ans2)
                            # Longer span
                            d = float(input("What is the depth of the longer span: "))
                            # midspan
                            d_l = round(pdict[p] - cop, 2)  # design load(pressure)
                            m1 = round(s4 * d_l * (lx ** 2), 2)
                            k = round(m1 / ((Fcu * b) * d ** 2), 2)
                            La = round(0.5 + math.sqrt(0.25 - (k / 0.9)), 2)
                            Z = La * d
                            Ans3 = m1 / (0.95 * Fy * Z)
                            Ans3 = ("The area of the mid span for the first panel under the longer span is: ", + Ans3)
                            print(Ans3)
                            # continuos edge
                            m2 = round(s3 * d_l * (lx ** 2), 2)
                            k = round(m2 / ((Fcu * b) * d ** 2), 2)
                            La = 0.5 + math.sqrt(0.25 - (k / 0.9))
                            Z = La * d
                            Ans4 = m2 / (0.95 * Fy * Z)
                            Ans4 = (
                            "The area of the continuous edge for the first panel under the longer span is: ", + Ans4)
                            print(Ans4)
                            Answers.update({"Ans1": Ans1, "Ans2": Ans2, "Ans3": Ans3, "Ans4": Ans4})
                        else:
                            print("the ratio ly to lx is not less than or equal to two, hence enter a valid number ")
                    else:
                        print('alright')
                        span = round(float(input("input the span over the edges  ")), 3)
                        coef = round(float(input("Input the bending moment coefficient ")), 3)
                        w4 = Pc
                        Fcu4 = float(input("Input the Fcu "))
                        b = float(input("Input the breadth "))
                        m4 = coef * (span ** 2) * w4
                        k4 = round(m4 / (Fcu4 * b * (d ** 2)), 2)
                        L4 = round(0.5 + math.sqrt(0.25 - (k4 / 0.9)), 2)
                        z4 = round(L4 * d, 2)
                        Fy = float(input("Input the Fyv value "))
                        As4 = m4 / (L4 * Fy * z4)
                        anss = (0.15 / 100) * b * d
                        print(anss)
                        anss = ("The area of the continuous edge for the first panel under the longer span is: ", +  anss)
                        print(anss)
                        print(i)
                        Answers.update({"ans": anss})
            else:
                print('Enter the input again ')
            for k in Answers.values():
                print(k)
            return Answers
    def beam(self):
        print('Program for the development of a Raft Beam ')
        print("For the first gridline, input the follwoing data in m ")
        #a = input("How many beams are in this structure? ")
        print("Input all values in metres")
        lx = float(input("What is the lx "))
        ly = float(input("What is the ly "))
        wab = float(input("What is the weight of the span on beam AB "))
        wbc = float(input("What is the weight of the span on beam BC "))
        wcd = float(input("What is the weight of the span on beam CD "))
        h = float(input("What is the thickness of the slab "))
        b = float(input("What is the breadth of the slab "))
        gk = h * b * 24
        bol = 1.4 * gk
        sl1 = (1 / 3) * lx * gk
        sl2 = (1 / 3) * ly * gk
        sab = wab + bol + sl1
        sbc = wbc + bol + sl2
        scd = wcd + bol + sl2
        print(sab, sbc, scd)
        # h = design_load + bol + sl
        # Support Moment
        # Using three monent equation
        # Span moment
        mma = float(input("What is the load on AB "))
        mmb = float(input("What is the load on BC "))
        mmc = float(input("What is the load on CD "))
        mmab = float(input("What isvthe distance between A and B "))
        mmbc = float(input("What is the distance between B and C "))
        mmcd = float(input("What is the distance between C and D "))
        #to solve for the upward moment
        a1 = 2 * (mmab + mmbc)
        b1 = mmbc
        c1 = ((mma*(mmab**3)/4) + (mmb*(mmbc**3)/4))
        b2 = 2 * (mmbc + mmcd)
        a2 = mmcd
        c2 = ((mmb*(mmbc**3)/4) + (mmc*(mmcd**3)/4))
        A = np.array([[a1,b1], [a2,b2]])
        B = np.array([c1,c2])
        D = np.linalg.inv(A)
        E = np.dot(D,B)
        print(E)
        print("Since Ma = Md = 0")
        Ma = 0
        Md = 0
        Mb = E[0]
        Mc = E[1]
        Mab = (1 / 8) * sab * (lx ** 2)
        Mbc = (1 / 8) * sbc * (ly ** 2)
        Mcd = (1 / 8) * scd * (ly ** 2)
        print(Mab, Mbc, Mcd)
        # due to superimposition
        mab = Mab - (Mb / 2)
        mbc = Mbc - ((Mb + Mc) / 2)
        mcd = Mcd - (Mc / 2)
        mmm = [mab, mbc, mcd]
        # design #Span A-B
        # span = len(mmm)
        for i in mmm:
            M = i #float(input('Input the moment of the span AB in mm '))
            b = float(input('Input the breadth in mm '))
            d = float(input('Input the height in mm '))
            Fcu = float(input('Input the Fcu '))
            k = M / (Fcu * b * (d ** 2))
            La = round(0.5 + math.sqrt(0.25 - (k / 0.9)), 2)
            if La <= 0.95:
                pass
            else:
                La = 0.95
            Z = La * d
            print(Z, "mm")
            Fy = float(input("What is your Fy "))
            As = M / (0.95 * Fy * Z)
            Minimum_steel = 0.15 * b * d
            print(As, "mm^2")
            print(Minimum_steel, "mm^2")
            #shear
            print("The Shear force on the respective beams are; ")
            Vab = ((sab * lx) / 2) + ((Ma - Mb) / lx)
            print(Vab, "KN")
            Vba = ((sab * lx) / 2) + ((Mb - Ma) / lx)
            print(Vba, "KN")
            Vbc = ((sbc * ly) / 2) + ((Mb - Mc) / ly)
            print(Vbc, "KN")
            Vcb = ((sbc * ly) / 2) + ((Mc - Mb) / ly)
            print(Vcb, "KN")
            Vcd = ((scd * ly) / 2) + ((Mc - Md) / ly)
            print(Vcd, "KN")
            Vdc = ((scd * ly) / 2) + ((Md - Mc) / ly)
            print(Vdc, "KN")
            #According to the shear stress above
            print("The shear stresses are; ")
            bd = 288000
            vab = Vab / (bd)
            print(vab, "N/mm^2")
            vba = Vba / (bd)
            print(vba, "N/mm^2")
            vbc = Vba / (bd)
            print(vbc, "N/mm^2")
            vcb = Vcb / (bd)
            print(vcb, "N/mm^2")
            vcd = Vcd / (bd)
            print(vcd, "N/mm^2")
            vdc = Vdc / (bd)
            print(vdc, "N/mm^2")
            #Grid lines
            Asv = float(input("what is the value of Asv "))
            Fyv = float(input("What is the value of Fyv "))
            br = float(input("What's the value of br "))
            Sva = (Asv * 0.95 * Fyv) / (0.4 * br)
            print(Sva)
            Svb = (Asv * 0.95 * Fyv) / (br * (Vba - Vbc))
            print(Svb)
            #Gridline B
            print("For the data on the Gridline B, i.e. Beam B.1-4, input the following values in mm or N/mm ")
            ww = float(input("What is the total load on the gridline "))
            lxx = float(input("what is the lx "))
            lyy = float(input("What is the ly "))
            kk = ly / lx
            loading = 0.5 * ww * lxx * (1 - (1 / (3 * (kk ** 2))))
            ww1 = loading * 2
            ww2 = ww
            ww3 = ww1
            #thus total load on beam including own weight
            www1 = ww3 + bol + loading
            www2 = ww2 + bol + gk
            #Moment Analysis
            mmmab = float(input("What is the load on AB "))
            mmmbc = float(input("What is the load on BC "))
            mmmcd = float(input("What is the load on CD "))
            mmmabb = float(input("What isvthe distance between A and B "))
            mmmbbc = float(input("What is the distance between B and C "))
            mmmcdd = float(input("What is the distance between C and D "))
            #to solve for the upward moment
            aa1 = 2 * (mmmabb + mmmbbc)
            bb1 = mmmbbc
            cc1 = ((mmmabb*(mmmabb**3)/4) + (mmmbbc*(mmmbbc**3)/4))
            bb2 = 2 * (mmmbbc + mmmcdd)
            aa2 = mmmcdd
            cc2 = ((mmmbc*(mmmbc**3)/4) + (mmmcd*(mmmcdd**3)/4))
            AA = np.array([[aa1,bb1], [aa2,bb2]])
            BB = np.array([cc1,cc2])
            DD = np.linalg.inv(A)
            EE = np.dot(D,B)
            print(EE)
            Mm1 = 0
            Mm4 = 0
            Mm2 = EE[0]
            Mm3 = EE[1]
            #Span moment
            msab = (1/8) * www1 * lyy
            msbc = (1 / 8) * www2 * lxx
            mscb = msab
            print(msab, msbc, mscb)
            Mssab = msab - (Mm2 / 2)
            Mssbc = msbc - ((Mm2 + Mm3) / 2)
            Msscb = mscb - (Mm2 / 2)
            print(Mssab, Mssbc, Msscb)
            #design
            kkk = ((Mssab * (10 ** 2)) / (Fcu * b * (h ** 2)))
            LLa = round(0.5 + (math.sqrt(0.25 - (kkk / 0.9))), 4)
            zzz = LLa * d
            print(zzz, LLa, kkk)
            ans = ((Mssab * (10 ** 2)) / (0.95 * Fy * zzz))
            print(ans, "mm2")
            kkb = ((Mm3 * (10 ** 2)) / (Fcu * b * (d ** 2)))
            LLb = 0.5 + math.sqrt(0.25 - (kkb / 0.9))
            zzb = LLb * d
            print(zzb, LLb, kkk)
            anb = ((Mm3 * (10 ** 2)) / (0.95 * Fy * zzb))
            print(anb, "mm2")
            #shear stress
            Vvab = ((mmmab*mmmabb) / 2) + ((Mm1 - Mm2)/mmmabb)
            Vvba = ((mmmab*mmmabb) / 2) + ((Mm2 - Mm1)/mmmabb)
            Vvbc = ((mmmbc*mmmbbc) / 2) + ((Mm2 - Mm3)/mmmbbc)
            Vvab = Vvdc
            Vvba = Vvcd
            Vvbc = Vvcb
            #Wab, Wba, Wbc =✓weight on beam
            #for maximum capacity
            ym = float(input("What's the value of ym "))
            vvc = ((0.79/ym) * (((100 * anss) ** (1/3)) / bd) * ((400 / d) ** (1/4)) * ((Fcu/25) ** (1/3)))
            print(vvc)
            print(vvc + 0.4)
            print(0.5 * vvc)
            #bbbbbbb #shear node
            vvab = ((Vvab * 1000) / bd)
            vvba = ((Vvba * 1000) / bd)
            vvbc = ((Vvbc * 1000) / bd)
            #note: if a span falls under the category (Vc + 0.4) < v < 0.8
            ssv = ((anss * 0.95 * Fyv) / b* (vvba - vvc))
            print(vvc)
            print("Note: When the spacing of links so calculated is below 125mm, it is advisable to double the sv and provide double the leg of the link reinforcement. as in the above, four legs of 10mm diameterbars at double sv, that's 200mm.")
  
a = raft_design()
print(a.panel())
print(a.beam())





