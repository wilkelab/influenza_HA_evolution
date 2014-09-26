import Bio.AlignIO, sys
import numpy as np
import scipy.stats as ss
from sets import Set
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

allowed = ['a','c','t','g']

gencode = { 'ATA':'I',
  'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T',
  'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S',
  'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L',
  'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H',
  'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R',
  'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A',
  'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E',
  'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S',
  'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L',
  'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C',
  'TGT':'C', 'TGA':'*', 'TGG':'W', }

max_sa = {'A':129 ,'R':274, 'N':195, 'D':193, 'C':167, 'E':223, 'Q':225, 'G':104, 'H':224, 'I':197, 'L':201, 'K':236, 'M':224, 'F':240, 'P':159, 'S':155, 'T':172, 'W':285, 'Y':263, 'V':174}

#Generate epitope dictionary
eps = {}
#for key in [122, 124, 126, 130, 131, 132, 133, 135, 137, 138, 140, 142, 143, 144, 145, 146, 150, 152, 168]:
#  eps[key] = 'A'
#for key in [128, 129, 155, 156, 157, 158, 159, 160, 163, 164, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198]:
#  eps[key] = 'B'
#for key in [44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312]:
#  eps[key] = 'C'
#for key in [96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248]:
#  eps[key] = 'D'
#for key in [57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265]:
#  eps[key] = 'E'
#for key in [49, 60, 74, 79, 90, 274, 151, 52, 277, 220, 134, 136, 153, 17, 199, 2, 3, 4, 31, 112, 205, 220, 271]:
#  eps[key] = 'O'
#for key in [98, 135, 136, 137, 138, 153, 155, 183, 190, 194, 195, 224, 225, 226]:
#  eps[key] = 'R'

for key in [18,  20,  34,  37,  38,  54,  66,  67,  98, 105, 106, 107, 108, 110, 115, 117, 119, 120, 121, 122, 124, 126, 127, 128, 129, 130, 131, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 155, 156, 157, 158, 159, 160, 173, 174, 175, 176, 177, 178, 179, 180, 181, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 222, 225, 226, 259, 260, 262, 276, 278, 289, 291, 318, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 347, 348, 349, 350, 363, 364, 366, 367, 368, 370, 371, 374, 375, 377, 378, 379, 381, 382, 385, 387, 388, 389, 482]:
    eps[key] = 'A'

def main():
    stats = importr('stats')
    
    args =  sys.argv
    file1 = args[1]
    file2 = args[2]
    
    frequencies1, rep_seq1 = calculate_frequencies(file1, False)
    frequencies2, rep_seq2 = calculate_frequencies(file2, False)
    
    freq_list, pvalue_list = calc_pvalues(frequencies1, frequencies2)
    
    pvalue_list = np.array(pvalue_list)
    pvals_only = pvalue_list[:, 1]
    
    #save_aa_freq_trajectories(freq_list, pvalue_list)
    
    p_adjust = stats.p_adjust(ro.FloatVector(pvals_only), method = 'BH')
    
    pvalue_list[:, 1] = p_adjust
    rsas = get_rsa('4FNK.pdb')
           
    rsa_select, num_selected, num_in_eps = output_positive_site(rsas, pvalue_list, rep_seq1)
    
    positive_sites_mean_rsa = np.mean(rsa_select)
    
    all_ep = []
    for value in set(eps.values()):
        rsa_ep = []
        for key in invert_dict(eps)[value]:
            if key in rsas.keys():
                rsa_ep.append(rsas[key][0])
                if value in ['A', 'B', 'C', 'D', 'E']:
                    all_ep.append(rsas[key][0])
    
    frac_real = 0
    if num_selected > 0:
        frac_real = float(num_in_eps) / num_selected

    values = np.array(rsas.values())[:,0]
    values = values.astype(float)
    
    output_other_crap(frac_real, all_ep, positive_sites_mean_rsa, rsa_select, values)
    
    #Calculate the chances of having the sites we found be in a random epitope with the same RSA distribution as the real epitopes
    pval_random = []
    for i in pvalue_list:
        if int(i[0]) in rsas.keys() and int(i[0]) in eps.keys():
            if rsas[int(i[0])] > positive_sites_mean_rsa:
                pval_random.append(i[1])
                #print(rsas[int(i[0])])
    
    ##Turn off bootstrapping to make the calculation go faster
    probability_random = calc_probability_epitope_random(10, rsas, pvalue_list, all_ep, frac_real)
    
    out_file = open('fraction_real_P.dat', 'a')
    out_file.write(str(frac_real) + '\t' + str(probability_random[0]) + '\t' + str(probability_random[1]) + '\n')
    out_file.close()

    print('Average p-value for high RSA sites: ' + str(probability_random))

def output_other_crap(frac_real, all_ep, positive_sites_mean_rsa, rsa_select, values):
    print("\nFraction positive sites in real epitopes: " + str(frac_real))
    print('Average RSA of all epitopes sites:       ' + str(np.mean(all_ep)))
    print('Average RSA of selected site:            ' + str(positive_sites_mean_rsa))
    print('Average RSA of all sites:                ' + str(np.mean(values)))
    print('RSA of selection versus RSA of epitopes: ' + str(ss.mannwhitneyu(rsa_select, all_ep)[1]))
    print('RSA of selection versus RSA of all:      ' + str(ss.mannwhitneyu(rsa_select, values)[1]))
    print('RSA of epitopes versus RSA of all:       ' + str(ss.mannwhitneyu(all_ep, values)[1]))
    
def output_positive_site(rsas, pvalue_list, rep_seq1):
    rsa_select = []
    num_in_eps = 0
    num_selected = 0 
    
    #print("\nEpitopes from the literature....\n")
    for i in pvalue_list:
        if i[1] < 0.001:
            out_file = open('positive_sites_P.dat', 'a')
            out_file.write(str(i[0]) + '\n')
            out_file.close()

        if i[1] < 0.001 and int(i[0]) in eps.keys() and int(i[0]) in rsas.keys():
            print(i[0], rep_seq1[int(i[0])+(16-1)], rsas[int(i[0])][1], eps[int(i[0])], rsas[int(i[0])][0], i[1])
            rsa_select.append(rsas[int(i[0])][0])
            num_selected += 1
            
            #See if the site is in an epitope
            if eps[int(i[0])] in ['A', 'B', 'C', 'D', 'E']:
                num_in_eps += 1
                
        elif i[1] < 0.001 and int(i[0]) in eps.keys():
            print(i[0], rep_seq1[int(i[0])+(16-1)], 'None', eps[int(i[0])], 'None', i[1])
            num_selected += 1
            
            #See if the site is in an epitope
            if eps[int(i[0])] in ['A', 'B', 'C', 'D', 'E']:
                num_in_eps += 1
                
        elif i[1] < 0.001 and int(i[0]) in rsas.keys():
            print(i[0], rep_seq1[int(i[0])+(16-1)], rsas[int(i[0])][0], 'None', rsas[int(i[0])][0], i[1])
            rsa_select.append(rsas[int(i[0])][0])
            num_selected += 1
        elif i[1] < 0.001:
            print(i[0], rep_seq1[int(i[0])+(16-1)], 'None', 'None', 'None', i[1])
            num_selected += 1
        #This line is just for QC to make sure the alignment is correct
        #elif int(i[0]) in eps.keys() and int(i[0]) in rsas.keys():
        #    #print(i[0], rep_seq1[int(i[0])+(16-1)], rsas[int(i[0])][1])
    return(rsa_select, num_selected, num_in_eps)
    
def calc_pvalues(frequencies1, frequencies2):
    stats = importr('stats')
    freq_list = []
    
    pvalue_list = []
    
    for obs in range(len(frequencies2)):
        temp_obs = np.array(frequencies2[obs].values())
        temp_exp = np.array(frequencies1[obs].values())
        
        freq_list.append(temp_exp/np.sum(temp_exp))

        if not np.array_equal(temp_obs/np.sum(temp_obs), temp_exp/np.sum(temp_exp)):
            equivalent = []
            for ex, ob in zip(temp_exp, temp_obs):
                if ex == 0 and ob == 0:
                    equivalent.append(False) 
                else: 
                    equivalent.append(True)
                    
            equivalent = np.array(equivalent)
             
            #Using scipy chi-squared function doesn't allow simulating p-values     
            #tmp_ss = ss.chi2_contingency([temp_obs[equivalent], temp_exp[equivalent]])
            #pvalue_list.append([obs - 16 + 1, tmp_ss[1]])
            
            #Using rpy chi-squared function does allow bootstrapping p-values
            tmp_matrix = ro.r.matrix(ro.FloatVector(temp_obs[equivalent]) +  ro.FloatVector(temp_exp[equivalent]), nrow=2)
            tmp_matrix = tmp_matrix.transpose()
            tmp_ss = stats.chisq_test(tmp_matrix, simulate_p_value=True, B=100000)[2][0]
            pvalue_list.append([obs - 16 + 1, tmp_ss])
        else:
            pvalue_list.append([obs - 16 + 1, 1])
            
    return(freq_list, pvalue_list)
            
def mean_confidence_interval(data, confidence=0.999):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), ss.sem(a)
    h = se * ss.t.ppf((1+confidence)/2., n-1)
    return m, h
       
def invert_dict(d):
    newdict = {}
    for k, v in d.iteritems():
        newdict.setdefault(v, []).append(k)
    return newdict
    
def calculate_frequencies(in_file, freqs):
    alignment = Bio.AlignIO.read(in_file, 'fasta')

    dat = []
    for sequence in alignment:
        temp_dat = []
        for nuc in range(len(sequence.seq)/3 - 1):  
            codon = sequence.seq[(nuc*3)] + sequence.seq[(nuc*3)+1] + sequence.seq[(nuc*3)+2]
            temp_dat.append(codon)
        dat.append(temp_dat)
    easy_dat = np.array(dat)

    aas = []
    rep_seq = []
    for count in range(len(easy_dat[0,:])):
        temp_aas = []
        #Don't try to understand this
        codons = [temp for temp in easy_dat[:, count] if temp[0] in allowed and temp[1] in allowed and temp[2] in allowed]
        for codon in range(len(codons)):
            temp_aas.append(gencode[codons[codon].upper()])
            if codon == 0:
                rep_seq.append(gencode[codons[codon].upper()])
        aas.append(np.array(temp_aas))

    frequencies = []

    for site in aas:
        counts = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}
        types = set(site)
        for aa in types:
            if freqs and aa in counts.keys():
              counts[aa] = float(list(site).count(aa))/len(site)
            elif aa in counts.keys():
              counts[aa] = float(list(site).count(aa))

        frequencies.append(counts)

    return(np.array(frequencies), rep_seq)

def get_rsa(pdb_file):
    import subprocess, os
    import dssp_parse as dd
    
    subprocess.call('mkdssp -i ' + pdb_file + ' -o tmp.dssp', shell=True)
    
    dd_ob = dd.DSSPData()
    dd_ob.parseDSSP('tmp.dssp')
    os.remove('tmp.dssp')
    
    accs = dd_ob.getACC()
    aas = dd_ob.getAAs()
    resnums = range(9,326) + range(330, 503)
    rsas = {}
    
    for aa, acc, resnum in zip(aas, accs, resnums):
        if aa in max_sa.keys():
            rsas[int(resnum)] = [float(acc)/max_sa[aa], aa]
    return(rsas)
    
def calc_prob_dist(values, all_ep, num_bins):
    data = np.array(all_ep)
    bins = np.arange(0., 1.0, float(1)/num_bins)
    
    bin_counts = np.zeros(num_bins)
    bin_counts += 1
    
    for i in data:
        for j in range(len(bins) - 1):
           if bins[j] <= i and i < bins[j + 1]:
               bin_counts[j] += 1
        if i >= 1:
            bin_counts[num_bins - 1] += 1
    
    bin_counts = bin_counts/np.sum(bin_counts)
    prob_count = np.zeros(num_bins)
    prob_count += 1
    
    for i in values:
        for j in range(len(bins) - 1):
            if bins[j] <= i[1] and i[1] < bins[j + 1]:
                prob_count[j] += 1
        if i[1] >= 1:
            prob_count[num_bins - 1] += 1
    
    probabilities = []
    for i in values:
        for j in range(len(bins) - 1):
            if bins[j] <= i[1] and i[1] < bins[j + 1]:
                probabilities.append(1 / prob_count[j] * bin_counts[j])
        if i[1] >= 1:
            probabilities.append(1 / prob_count[num_bins - 1] * bin_counts[num_bins - 1])
    
    probabilities = probabilities/np.sum(probabilities)
    random_sites = np.random.choice(range(len(values)), len(data), p = probabilities)
    
    return(np.array(values)[random_sites])
    
def calc_probability_epitope_random(num_bins, rsas, pvalue_list, all_ep, frac_real):

    ##################################################################################################################
    #Start to figure out the probability that a site will be in an epitope given a random sample of similar RSA sites#
    ##################################################################################################################
    
    values = np.array(rsas.values())[:,0]
    values = values.astype(float)
    sites = np.array(rsas.keys())
    recombined = [[site, value] for site, value in zip(sites, values)]
    
    probability_random = []
    
    
    for i in range(1000):
        rsa_prob_dist = calc_prob_dist(recombined, all_ep, num_bins)
    
        test_eps = np.array(rsa_prob_dist)[:,0]
    
        for key in eps.keys():
            if key in test_eps:
                eps[key] = 'Ep'
            else:
                eps[key] = 'NE'
            
        num_in_eps = 0
        num_selected = 0 
        rsa_select = []
    
        #print("\nRandom epitopes based on RSA....\n")
        for i in pvalue_list:
            if i[1] < 0.01 and int(i[0]) in eps.keys() and int(i[0]) in rsas.keys():
                #print(i[0], rep_seq1[int(i[0])+(16-1)], rsas[int(i[0])][1], eps[int(i[0])], rsas[int(i[0])][0], i[1])
                rsa_select.append(rsas[int(i[0])][0])

                num_selected += 1
                #See if the site is in an epitope
                if eps[int(i[0])] == 'Ep':
                    num_in_eps += 1
                
            elif i[1] < 0.01 and int(i[0]) in eps.keys():
                #print(i[0], rep_seq1[int(i[0])+(16-1)], 'None', eps[int(i[0])], 'None', i[1])
                num_selected += 1
            
                #See if the site is in an epitope
                if eps[int(i[0])] == 'Ep':
                    num_in_eps += 1
                
            elif i[1] < 0.01 and int(i[0]) in rsas.keys():
                #print(i[0], rep_seq1[int(i[0])+(16-1)], rsas[int(i[0])][0], 'None', rsas[int(i[0])][0], i[1])
                rsa_select.append(rsas[int(i[0])][0])
                num_selected += 1
            elif i[1] < 0.01:
                #print(i[0], rep_seq1[int(i[0])+(16-1)], 'None', 'None', 'None', i[1])
                num_selected += 1    
            
        frac_random = float(num_in_eps) / num_selected
        probability_random.append(frac_random)
    
    probability_random = np.array(probability_random)
    #print(np.sum([probability_random > frac_real]))
    prob = np.sum([probability_random > frac_real]) / float(len(probability_random))    
    #print("\nFraction positive sites in random epitopes: " + str(prob) + "\n")
    
    #print(mean_confidence_interval(probability_random))
    
    outfile = open('prob_data_P.dat', 'a')
    line = ''
    for i in probability_random:
        line += ',' + str(i)
    outfile.write(line[1:len(line)] + '\n')
    outfile.close()
    
    return(mean_confidence_interval(probability_random))

def save_aa_freq_trajectories(freq_list, pvalue_list):
    #Save the amino acid frequency trajectories
    for i in range(len(freq_list)):
        output = open('trajectories/traj_' + str(int(pvalue_list[i][0])) + '.txt', 'a')
        line = ''
        for j in freq_list[i]:
            line += str(j) + ','
        output.write(line[0:len(line)-1] + '\n')
        output.close()
        
main()
