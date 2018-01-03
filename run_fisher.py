import sys
sys.path.insert(0, 'source')
import fisher_which

""" Doing fisher for cmb or 21cm """
#cmb_or_21 = '21'
cmb_or_21 = 'cmb'

""" Fixing Yp with BBN or not """
Yp_BBN = True
#Yp_BBN = False

fisher_which.run (cmb_or_21, Yp_BBN)
