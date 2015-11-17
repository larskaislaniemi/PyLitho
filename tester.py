import Crust
import PrplxWrap

components = ["NA2O", "MGO", "AL2O3", "SIO2", "K2O", "CAO", "TIO2", "MNO", "FEO", "H2O"]
inimasses = [0.0327, 0.0248, 0.1540, 0.6662, 0.0280, 0.0359, 0.0064, 0.0010, 0.0504, 0.0100]

inimasses = [2800.0*mass for mass in inimasses]

perplex = PrplxWrap.perplex("crustmod", callerCompnames=components)
cont = Crust.Container(compos=components, masses=inimasses, perplex=perplex)

cont.P = 13700.
cont.T = 950.

f = open('out.csv', 'a')

f.write('"P","T","dT","melt"')
f.write("\n")
cont.updatePerplex(perplex=perplex)
for i in range(100):
    prevT = cont.T
    cont.move_isentropic(-10, perplex=perplex)
    if "Melt" in cont.perplexResult['NAMEPHASES']:
        meltwt = cont.perplexResult['WTPHASES'][cont.perplexResult['NAMEPHASES'].index("Melt")]
    else:
        meltwt = 0.0
    print str(cont.P) + "," + str(cont.T) + "," + str(cont.T-prevT) + "," + str(meltwt)
    f.write(str(cont.P) + "," + str(cont.T) + "," + str(cont.T-prevT) + "," + str(meltwt))
    f.write("\n")


f.close()
