import numpy as np

amplifier_gains = np.array([26,23,15]) #Amplifier gains in dBs
amplifier_NFs = np.array([5.5,5.5,6]) #Amplifier noise figures in dBs
amplifier_sat_op_powers = np.array([20,17,20]) #Amplifier saturated ouput power in dBms

transceiver_mode_min_req_SNRs = {"A":21.0,"B":13.8,"C":8.0,"D":2.6} #Transceiver minimum required SNRs in dBs based on mode chosen
transceiver_mode_min_req_opt_powers = {"A":-24.8,"B":-26.0,"C":-27.8,"D":-30.8} #Transceiver minimum required output power in dBms based on mode chosen
transceiver_net_data_rates = {"A":800,"B":600,"C":400,"D":200} #Transceiver net data rates in Gb/s based on mode chosen

traffic_vals_edges_ring_percent =  {"Lon-B":28,"B-M":18,"M-G":3,"G-Le":8,"Le-Lon":19} #Based on hand drawn diagram
traffic_vals_edges_meshA_percent = {"Lon-B":12,"B-M":6,"Lon-M":19,"M-G":3,"M-Le":4,"G-Le":1,"Le-Lon":8}
traffic_vals_edges_meshB_percent = {"Lon-B":12,"B-M":4,"Lon-M":12,"Lon-G":7,"M-G":3,"B-Le":2,"M-Le":2,"G-Le":1,"Le-Lon":8}

#0 for bypass choice means the automated choice is used (it is not bypassed) - less relevant for current setup
topology_values = {"ring":[160,traffic_vals_edges_ring_percent,["Lon-B","B-M","M-G","G-Le","Le-Lon"],3],"meshA":[236,traffic_vals_edges_meshA_percent,["Lon-B","B-M","Lon-M","M-G","M-Le","G-Le","Le-Lon"],3],"meshB":[375,traffic_vals_edges_meshB_percent,["Lon-B","B-M","Lon-M","Lon-G","M-G","B-Le","M-Le","G-Le","Le-Lon"],3]}

edge_proper_names = {"Lon-B":"London-Birmingham","B-M":"Birmingham-Manchester","Lon-M":"London-Manchester","Lon-G":"London-Glasgow","M-G":"Manchester-Glasgow","B-Le":"Birmingham-Leeds","M-Le":"Manchester-Leeds","G-Le":"Glasgow-Leeds","Le-Lon":"Leeds-London"}

R_s = 95 #transceiver_symbol_rate in GBd
R_s_SI = R_s * 10**9
alpha_dB_per_km = 0.17 #global as fixed for a fibre
alpha_per_km = np.log(10**(0.1*alpha_dB_per_km)) #global as fixed for a fibre
alpha_per_m = alpha_per_km*(10**(-3))
D_fibre = 18 #fixed worst case dispersion in ps/(nm km), can change value or use range() to form array to test range 
#of values
beta_2 = abs(-D_fibre/0.78) #beta_2 in (ps)^2 / km which is the second derivative of the with respect to the angular 
#frequency, again worst case => can change to form range
beta_2_SI = beta_2 * 10**(-27)
BW_mod_total = 4.5 * 10**12 #Total modulated bandwidth as found by 45 wavelengths spaced at 100 GHz

h = 6.63 * 10**(-34)
lambda_0 = 1550 * 10**(-9) #fixed C band centre wavelength in m (around which WDM occurs)
nu_0 = 193.1 * 10**12#(2.998*(10**8))/lambda_0

n2 = 2.6 * 10 **(-20) #constant for Kerr effect
k_0 = (2*np.pi)/lambda_0 #wavenumber 
A_eff = 85 * 10**(-12) #Effective area of fibre in m^2
gamma_kerr = (n2*k_0)/(A_eff)

#No docstrings for the first four functions as self-evident from name

def dB_to_lin(dB_val): #Convert from dB to linear value for power
    return 10**(0.1*dB_val)

def lin_to_dB(lin_val): #Convert from linear to dB value for power
    return 10*np.log(lin_val)/np.log(10)

def dBm_to_lin(dBm_val): #Convert from linear to dBm value for power
    return 10**(0.1*dBm_val) * 10**(-3) 

def lin_to_dBm(lin_val): #Convert from dBm to linear value for power
    return 10*np.log(lin_val/(10**(-3)))/np.log(10)

#The initial part of this function remians relevant while the later part is redundant as the choice to minimise N_ASE 
#(i.e. EDFA3) should always be chosen for SNR and thus, the bypass option is usually used
def amplifier_choice(choices,L,secondary_link,amp_saturated,bypass_choice):
    """
    Input parameters:
    choices: EDFA numbers available
    L: span length in km 
    secondary_link: boolean value to check if connected to transmitter node - if so, False (account for ROADM during
    amplifier choice) else True
    amp_saturated: boolean to check if EDFA is operating under gain saturation
    bypass_choice: bypass second part (NF and saturation output power) of EDFA choice process so that same EDFA 
    chosen in every case as long as it satisfies first part (gain check)
    
    Returns:
    EDFA_no: EDFA number chosen (1,2,3 or 0 if no appropriate options - will prompt user to choose different 
    distances in a chronologically later function)
    """
    choices_dummy = choices.copy() #Going to edit choices to eliminate them => work on copy
    for ind,choice in enumerate(choices_dummy): #Check to see if amplifier cannot compensate for loss at given value for L
        corrected_gain = amplifier_gains[choice-1] - 3 #Assume 3 dB loss when choosing amplifier
        if (secondary_link == False & amp_saturated == False):
            if (alpha_dB_per_km*L + 4 > amplifier_gains[choice-1]):
                choices_dummy[ind] = 0
        elif (secondary_link == False & amp_saturated == True):
            if (alpha_dB_per_km*L + 4 > corrected_gain):
                choices_dummy[ind] = 0
        elif (secondary_link == True & amp_saturated == False):
            if (alpha_dB_per_km*L > amplifier_gains[choice-1]):
                choices_dummy[ind] = 0
        else:
            if (alpha_dB_per_km*L > corrected_gain):
                choices_dummy[ind] = 0
    choices_dummy_remaining_after_gain_check = [i for i in choices_dummy if i!= 0] #Remove choices_dummy that cannot meet gain requirement

    if (bypass_choice != 0 and bypass_choice in choices_dummy_remaining_after_gain_check): #Ensure that bypass only occurs if gain relation satisfied
        return bypass_choice

    if len(choices_dummy_remaining_after_gain_check)==1:
        EDFA_no = choices_dummy_remaining_after_gain_check[0]
        return EDFA_no
    elif len(choices_dummy_remaining_after_gain_check)==0:
        return 0
    #print(choices_dummy_remaining_after_gain_check) #for debugging

    #Set up so that noise figure takes precedence over saturation output power because will check for gain 
    #saturation to eliminate chocies anyway

    amplifier_NFs_avaiable = [NF for ind,NF in enumerate(amplifier_NFs) if ind+1 in choices_dummy]
    min_NF = np.amin(amplifier_NFs_avaiable)
    best_NF_choices_dummy = [ind+1 for ind,NF in enumerate(amplifier_NFs) if NF==min_NF]
    choices_dummy_remaining_after_NF_check = [choice for choice in choices_dummy_remaining_after_gain_check if choice in best_NF_choices_dummy]#Redundant if EDFA1, EDFA2 choice but this allows for the general case where EDFA 1 and 2 would not 
    #have the same NF
    #print(choices_dummy_remaining_after_NF_check) #for debugging
    #EDFA 3 usually eliminated because of worst gain and worst noise figure
    if len(choices_dummy_remaining_after_NF_check)==1:
        EDFA_no = choices_dummy_remaining_after_NF_check[0]
        return EDFA_no
    
    #Deciding factor is saturation output power if above does not return
    amplifier_sat_op_powers_available = [P for ind,P in enumerate(amplifier_sat_op_powers) if ind+1 in choices_dummy]
    max_P_out_sat = np.amax(amplifier_sat_op_powers_available)
    best_gain_sat_choices_dummy = [ind+1 for ind,P in enumerate(amplifier_sat_op_powers) if P==max_P_out_sat]
    choices_dummy_remaining_now = [choice for choice in choices_dummy_remaining_after_NF_check if choice in best_gain_sat_choices_dummy]
    #print(choices_dummy_remaining_now) #for debugging

    EDFA_no = choices_dummy_remaining_now[0] #if equal, then just choose either, set to first index
    return EDFA_no

def SNR_link_initlal_calc(link_name,EDFA_choices,L,secondary_link,bypass_automated_EDFA_choice):
    """
    Input parameters:
    link_name: name of link (usually includes city names) for the sake of EDFA number printed message used in past 
    iterations for debugging
    EDFA_choices: EDFA number
    L: span length in km 
    secondary_link: boolean value to check if connected to transmitter node - if so, False (account for ROADM during
    amplifier choice) else True
    bypass_automated_EDFA_choice: EDFA number that has been chosen instead of automated choice
    
    Returns:
    SNR_dBs: signal to noise ratio for the link in dBs
    G_opt: optimum power spectral density for the link in SI units
    N_ASE: amplified spontaneous emission noise from amplifier in SI units
    amp_choice: amplifier choice based on earlier function 
    """
    if (L==0):
        return 0,0,0,0 #So that additional links can be left unused without affecting operation
    L_eff = (1-np.exp(-alpha_per_km*L))/(alpha_per_km) * 10**3 #Effective length in m 
    C_NLI = (8*(gamma_kerr**2)*(L_eff**2)*alpha_per_m)/(27*np.pi*beta_2_SI)*np.log((beta_2_SI*(np.pi**2)*(BW_mod_total**2))/(alpha_per_m)) #Input everything into C_NLI eqn in SI units
    amp_choice = amplifier_choice(EDFA_choices,L,secondary_link,False,bypass_automated_EDFA_choice) #Choose amplfiier based on link length/span 
    amp_gain = amplifier_gains[amp_choice-1]
    amp_NF = amplifier_NFs[amp_choice-1] #Amplifier noise figure in dB
    N_ASE = dB_to_lin(amp_NF)*h*nu_0*(dB_to_lin(amp_gain)-1) #N_ASE equation based on initial amplifier choice 
    #N_ASE fixed by amplifier choice (usually EDFA 3)
    G_opt = (N_ASE/(2*C_NLI))**(1/3) #Expression for optimum power spectral density which can be used to find optimum launch power in SI units
    SNR_lin = 45*G_opt/(N_ASE*1.5) #Find SNR assuming 45 channels/wavelngths being used
    SNR_dBs = lin_to_dB(SNR_lin)
    #print("EDFA {} chosen for {} link".format(amp_choice,link_name)) # for debugging
    return SNR_dBs,G_opt,N_ASE,amp_choice

def operation_under_gain_saturation(link_name,EDFA_choices,L,secondary_link,WDM_lambdas_being_used,P_out_one_lambda,bypass_automated_EDFA_choice):
    """
    Input parameters:
    link_name: name of link (usually includes city names) for the sake of EDFA number printed message being 
    EDFA_choices: EDFA number
    L: span length in km 
    secondary_link: boolean value to check if connected to transmitter node - if so, False (account for ROADM during
    amplifier choice) else True
    WDM_lambdas_being_used: number of channels/wavelengths being used by link for wavelength division multiplexing
    P_out_one_lambda: optical power transmitted when using WWDM channel/wavelength used if actual amplifier gain used
    rather than -3dB reduced value in linear units (W)
    
    Returns:
    SNR_dBs: signal to noise ratio for the link in dBs
    break_condition: boolean to see if break condition satisfied
    """
    if (L==0):
        return 0 #So that additional links can be left unused without affecting operation
    break_condition = False
    L_eff = (1-np.exp(-alpha_per_km*L))/(alpha_per_km) * 10**3 #Effective length in m 
    C_NLI = (8*(gamma_kerr**2)*(L_eff**2)*alpha_per_m)/(27*np.pi*beta_2_SI)*np.log((beta_2_SI*(np.pi**2)*(BW_mod_total**2))/(alpha_per_m)) #Input everything into C_NLI eqn in SI units
    amp_choice = amplifier_choice(EDFA_choices,L,secondary_link,False,bypass_automated_EDFA_choice) #Choose amplfiier based on link length/span 
    if (amp_choice == 0):
        break_condition == True
    amp_NF = amplifier_NFs[amp_choice-1] #Amplifier noise figure in dB
    amp_P_sat = dBm_to_lin(amplifier_sat_op_powers[amp_choice-1]) / np.log(2)
    P_out = P_out_one_lambda*WDM_lambdas_being_used
    corrected_amp_gain =  amplifier_gains[amp_choice-1]-3 
    corrected_amp_gain_real =  lin_to_dB(dB_to_lin(amplifier_gains[amp_choice-1]) * np.exp(-(P_out/amp_P_sat)))
    #print("Gains comapred: {}, {}".format(corrected_amp_gain,corrected_amp_gain_real))

    #Change corrected_amp_gain to corrected_amp_gain_real to avoid using 3dB reduction in gain approximation 
    #(assumption)

    N_ASE = dB_to_lin(amp_NF)*h*nu_0*(dB_to_lin(corrected_amp_gain)-1) #N_ASE equation based on initial amplifier choice
    G_opt = (N_ASE/(2*C_NLI))**(1/3) #Expression for optimum power spectral density which can be used to find optimum launch power
    SNR_lin = (WDM_lambdas_being_used*G_opt)/(N_ASE*1.5) #Find SNR 
    SNR_dBs = lin_to_dB(SNR_lin)
    print("EDFA {} operating under gain saturation chosen for {} link".format(amp_choice,link_name) +  "\n" +  "\n")
    return SNR_dBs,break_condition

def dealing_with_gain_sat(amp_choice,EDFA_choices,G_opt,N_ASE,WDM_lambdas_being_used,L,secondary_link,bypass_automated_EDFA_choice):
    """
    Input parameters:
    amp_choice: EDFA (amplifier) number
    EDFA_choices: EDFA number
    G_opt: Optimum power spectral density from original calculation (without gain saturation)
    N_ASE: amplified spontaneous noise from EDFA from original calculation (without gain saturation)
    WDM_lambdas_being_used: number of channels/wavelengths being used by link for wavelength division multiplexing
    L: span length in km 
    secondary_link: boolean value to check if connected to transmitter node - if so, False (account for ROADM during
    amplifier choice) else True
    
    Returns:
    G_opt: New value for optimum power spectral density 
    N_ASE: New value for amplified spontaneous noise from EDFA 
    amp_choice: mew EDFA number chosen when gain saturation is taken into account
    operate_under_gain_sat: boolean that will trigger operation with corrected gain if no EDFA that can operate 
    without gain saturation
    """
    operate_under_gain_sat = False
    EDFA_choices_dummy = EDFA_choices.copy() #Going to edit choices to eliminate them => work on copy
    L_eff = (1-np.exp(-alpha_per_km*L))/(alpha_per_km) * 10**3 #Effective length in m 
    C_NLI = (8*(gamma_kerr**2)*(L_eff**2)*alpha_per_m)/(27*np.pi*beta_2_SI)*np.log((beta_2_SI*(np.pi**2)*(BW_mod_total**2))/(alpha_per_m)) #Input everything into C_NLI eqn in SI units
    #print("DEBUG: {},{},{}".format(C_NLI,N_ASE,G_opt))
    amp_sat_op_power = amplifier_sat_op_powers[amp_choice-1] #Amplfier saturated output power in dBm
    while (G_opt*R_s_SI*WDM_lambdas_being_used >= dBm_to_lin(amp_sat_op_power)): #Check for gain saturation case
        print("Optical power P_opt: {}".format(lin_to_dBm(G_opt*R_s_SI*WDM_lambdas_being_used))) #for debugging -3dB assumption error
        info_string = "If attempting to operate without gain saturation, the current optimal choice, EDFA " + str(amp_choice) + " as per gain and NF, was eliminated due to gain saturation " +  "(may end up using this choice under gain saturation if no such choices)"
        print(info_string + "\n") #Print statement so that situation known
        EDFA_choices_dummy.remove(amp_choice) #Remove choice where gain saturation has occured to a significant extent 
        #(gain reduced by at least factor of 2)
        if len(EDFA_choices_dummy) == 0: #If no choices in spite of gain saturation
            operate_under_gain_sat = True
            break
        updated_choices = EDFA_choices_dummy
        new_amp_choice = amplifier_choice(updated_choices,L,secondary_link,True,bypass_automated_EDFA_choice) #Choose new amplifier which will meet needs without gain saturation
        if ((bypass_automated_EDFA_choice != 0) or (new_amp_choice == 0)): #Check for bypass case or no suitable amplifier without gain saturation case
            operate_under_gain_sat = True
            break
        amp_sat_op_power = amplifier_sat_op_powers[new_amp_choice-1] #Amplfier saturated output power in dBm
        amp_choice = new_amp_choice #Update amplifier choice
        amp_NF = amplifier_NFs[amp_choice-1] #Amplifier noise figure in dB
        #Below calculation not exact but good approx for EDFA1 and EDFA3 after examining presumed P_tx_opt as 
        #gain approxiamtely halves (-3dB)
        #Not for EDFA2 as gain even lower but gets eliminated anyway due to lower saturated output power
        amp_gain =  amplifier_gains[amp_choice-1]
        N_ASE = dB_to_lin(amp_NF)*h*nu_0*(dB_to_lin(amp_gain)-1) #N_ASE equation based on initial amplifier choice
        G_opt = (N_ASE/(2*C_NLI))**(1/3) #Expression for optimum power spectral density which can be used to find optimum launch power
        
        #P_tx_opt_single_lambda = G_opt * R_s_SI 
        #P_tx_opt = P_tx_opt_single_lambda*WDM_lambdas_being_used #Can be used to assign condition output to variable
    else:
        print("Optical power P_opt: {}".format(lin_to_dBm(G_opt*R_s_SI*WDM_lambdas_being_used)))#for debugging -3dB assumption error & transceiver optical power tables
        return G_opt,N_ASE,amp_choice,operate_under_gain_sat
    #Break condition response
    return 0,0,amp_choice,operate_under_gain_sat

def SNR_gain_sat_accounted(amp_choice,link_name,EDFA_choices,G_opt,N_ASE,WDM_lambdas_being_used,L,secondary_link,P_out,bypass_automated_EDFA_choice):
    """
    Input parameters:
    amp_choice: EDFA (amplifier) number
    link_name: name of link (usually includes city names) for the sake of EDFA number printed message being 
    EDFA_choices: EDFA number
    G_opt: Optimum power spectral density from original calculation (without gain saturation)
    N_ASE: amplified spontaneous noise from EDFA 
    WDM_lambdas_being_used: number of channels/wavelengths being used by link for wavelength division multiplexing
    L: span length in km 
    secondary_link: boolean value to check if connected to transmitter node - if so, False (account for ROADM during
    amplifier choice) else True
    P_out: optical power transmitted when using WWDM channel/wavelength used if actual amplifier gain used 
    rather than -3dB reduced value in linear units (W)
    bypass_automated_EDFA_choice: EDFA number that has been chosen instead of automated choice
    
    Returns:
    SNR_dBs: signal to noise ratio for the link in dBs
    """
    if (L==0):
        return 0 #So that additional links can be left unused without affecting operation
    
    break_condition = False

    G_opt,N_ASE,amp_choice,operate_under_gain_sat = dealing_with_gain_sat(amp_choice,EDFA_choices,G_opt,N_ASE,WDM_lambdas_being_used,L,secondary_link,bypass_automated_EDFA_choice)
    
    if operate_under_gain_sat == True:
        SNR_dBs_gain_sat, break_condition = operation_under_gain_saturation(link_name,EDFA_choices,L,secondary_link,WDM_lambdas_being_used,P_out,bypass_automated_EDFA_choice)
        if break_condition == True:
            print("Error has occured as no suitable optical amplifier. Change link lengths" + "\n")
            return 0
        else:
            return SNR_dBs_gain_sat
    else:
        SNR_lin = (WDM_lambdas_being_used*G_opt)/(N_ASE*1.5) #Find SNR 
        SNR_dBs = lin_to_dB(SNR_lin)
        print("EDFA {} chosen for {} link".format(amp_choice,link_name) + "\n" +  "\n")
        return SNR_dBs

def NSRs_addition(SNRs):
    """
    Input parameters:
    SNRs: array of SNRs to be combined

    Returns:
    NSRs_sum: linear value for sum of noise to signal ratios 
    """
    #if (np.count_nonzero(SNRs) == 1): #Would have been required for earlier code (not for current one), useful 
    #function for the future
        #return dB_to_lin(SNRs[0])**(-1) 
    SNRs_lin = [dB_to_lin(i) for i in SNRs if i != 0]
    NSRs = np.divide(([1]*len(SNRs_lin)),(SNRs_lin))
    NSRs_sum = np.sum(NSRs)
    return NSRs_sum

def transceiver_mode_choice(P_tx_opt,mode_choices,SNR):
    """
    Input parameters:
    P_tx_opt: optical power from transmitter in linear units (W)
    mode_choices: "A","B","C","D" to dictacte transceiver operation parameters
    SNR: overall SNR at RX (transceiver)

    Returns:
    transceiver_mode: the chosen mode mainly based on the SNR at RX
    """
    for mode in mode_choices:
        if (P_tx_opt < dBm_to_lin(transceiver_mode_min_req_opt_powers[mode])): #As discussed in the Q&A session, this is just to account for extreme cases 
            #(1 channel with super long link) and will not play a role usually
            pass
        else:
            if (SNR < transceiver_mode_min_req_SNRs[mode]):
                pass
            else:
                transceiver_mode = mode
                #print(transceiver_mode)
                break
    return transceiver_mode

def transceiver_operation(P,SNR):
    """
    Input parameters:
    P: output optical power
    SNR: overall SNR at RX (transceiver)

    Returns:
    transceiver_mode: the chosen mode mainly based on the SNR at RX
    transceiver_net_data_rate: the data rate in Gb/s at which 
    """
    transceiver_mode = transceiver_mode_choice(P,list(transceiver_mode_min_req_SNRs.keys()),SNR) #Choose transceiver mode
    #transceiver_mode_min_req_SNR = transceiver_mode_min_req_SNRs[transceiver_mode] #Chosen mode's 
    #minimum required SNR if needed
    #transceiver_mode_min_req_opt_power = dBm_to_lin(transceiver_mode_min_req_opt_powers[transceiver_mode]) #Chosen 
    #mode's minimum required optical power if needed
    transceiver_net_data_rate = transceiver_net_data_rates[transceiver_mode] #Chosen mode's data rate
    return transceiver_mode,transceiver_net_data_rate

def capacity_calcs(total_transceiver_pairs,transceiver_net_data_rate,traffic_vals_dict,traffic_vals_key):
    """
    Input parameters:
    total_trasnceiver_pairs: as implied by name
    transceiver_net_data_rate: Data rate in Gb/s
    traffic_vals_dict: traffic values on a given edge in the chosen topology in % 
    traffic_vals_key: traffic in % for link speicifed by key name in above dictionary

    Returns:
    channels: number of WDM channels/wavelengths being used 
    actual_link_capacity: maximum capacity in link in Gb/s
    """
    channels = np.rint((total_transceiver_pairs * traffic_vals_dict[traffic_vals_key])/100)
    actual_link_capacity = transceiver_net_data_rate * channels
    return channels,actual_link_capacity

class City_Link:
    def __init__(self, link_name, EDFA_numbers, SNR, G_opt, N_ASE, amp_choice):
        self.link_name = link_name
        self.EDFA_numbers = EDFA_numbers
        self.SNR = SNR
        self.G_opt = G_opt
        self.N_ASE = N_ASE
        self.amp_choice = amp_choice

    def SNR_initial(self,EDFA_numbers,link_length,secondary_link,bypass_automated_EDFA_choice):
        self.SNR,self.G_opt,self.N_ASE,self.amp_choice = SNR_link_initlal_calc(self.link_name,EDFA_numbers,link_length,secondary_link,bypass_automated_EDFA_choice)
        self.P_opt_one_lambda = self.G_opt * R_s_SI


class City_Edge:
    def __init__(self, edge_code_name, edge_proper_name, edge_P_opt, initial_SNR, actual_SNR, initial_transceiver_mode, actual_transceiver_mode, WDM_lambdas_being_used,initial_capacity,actual_capacity):
        self.edge_code_name = edge_code_name
        self.edge_proper_name = edge_proper_name
        self.edge_P_opt = edge_P_opt
        self.initial_SNR = initial_SNR
        self.actual_SNR = actual_SNR
        self.initial_transceiver_mode = initial_transceiver_mode
        self.actual_transceiver_mode = actual_transceiver_mode
        self.WDM_lambdas_being_used = WDM_lambdas_being_used
        self.initial_capacity = initial_capacity
        self.actual_capacity = actual_capacity


max_number_of_links_in_edge = 15 #This ensures that even the longest longest link from London to Glasgow will 

def network_edge_parameters(EDFA_numbers,dists,topology):
    """
    Input parameters:
    EDFA_numbers: EDFA choices available (to make this setup more general in case more EDFAs are available)
    dists: two dimensional array of link distances i.e. distances between amplifier locations (and cities 
    if a link is connected to a city) within an edge 
    topology: type of topology being considered (can be ring, mesh type A or mesh type B)

    Returns:
    city_link_info_dict: dictionary consisting of multiple dictionaries 
    """
    edges = list(topology_values[topology][1].keys())

    bypass_automated_EDFA_choice = topology_values[topology][3]

    edge_objects = []

    for ind,edge in enumerate(edges): 

        if (len(dists[ind]) != 15): #To make this more convenient for the user
            dists[ind] = dists[ind] + [0]*(15-len(dists[ind]))

        y = City_Edge(edge,edge_proper_names[edge],0,0,0,"A","A",0,0,0)

        edge_link_objects = []
        edge_P_opt_one_lambda_max = 0 #replace with maximum as this deifnes the transceiver optical power

        for i in range(max_number_of_links_in_edge): #Make the initial SNR calculation for a link and store link objects in a list
            x = City_Link("{} {}".format(y.edge_proper_name,i+1), EDFA_numbers,0,0,0,0)
            x.SNR_initial(EDFA_numbers,dists[ind][i],False if (i == 0) else True,bypass_automated_EDFA_choice)
            edge_P_opt_one_lambda_max = max(edge_P_opt_one_lambda_max,x.P_opt_one_lambda)
            edge_link_objects.append(x)
        
        #Set edge object attributes to form final output dictionary parameters without accounting for gain saturation
        y.edge_P_opt = edge_P_opt_one_lambda_max
        edge_SNRs_list = [link.SNR for link in edge_link_objects]
        y.initial_SNR = lin_to_dB((NSRs_addition(edge_SNRs_list))**(-1))
        y.initial_transceiver_mode, edge_data_rate = transceiver_operation(y.edge_P_opt,y.initial_SNR)
        y.WDM_lambdas_being_used, y.initial_capacity = capacity_calcs(topology_values[topology][0],edge_data_rate,topology_values[topology][1],edge)

        
        edge_SNRs_actual_list = []

        for j in range(max_number_of_links_in_edge): #Recalculate SNR accounting for optical power effect on gain (i.e. gain saturation)
            link_SNR_actual = SNR_gain_sat_accounted(edge_link_objects[j].amp_choice,"{} {}".format(y.edge_proper_name,j+1),EDFA_numbers,edge_link_objects[j].G_opt,edge_link_objects[j].N_ASE,y.WDM_lambdas_being_used,dists[ind][j],False if (j == 0) else True,y.edge_P_opt,bypass_automated_EDFA_choice)
            edge_SNRs_actual_list.append(link_SNR_actual)

        #Set edge object attributes to form final output dictionary parameters accounting for gain saturation
        y.actual_SNR = lin_to_dB((NSRs_addition(edge_SNRs_actual_list))**(-1))
        y.actual_transceiver_mode,edge_data_rate_actual = transceiver_operation(y.edge_P_opt,y.actual_SNR)
        y.actual_capacity = capacity_calcs(topology_values[topology][0],edge_data_rate_actual,topology_values[topology][1],edge)[1]

        edge_objects.append(y)

    city_link_info_dict= {}

    dicts_list = []

    #Form list of dictionaries for each edge
    for edge_object in edge_objects: 
        dicts_list.append({"initial SNR":edge_object.initial_SNR,"actual SNR":edge_object.actual_SNR,"transceiver optical power":lin_to_dBm(edge_object.edge_P_opt * edge_object.WDM_lambdas_being_used),"initial transceiver mode":edge_object.initial_transceiver_mode,"actual transceiver mode":edge_object.actual_transceiver_mode,"number of WDM channels/wavelengths":edge_object.WDM_lambdas_being_used,"initial capacity":edge_object.initial_capacity,"actual capacity":edge_object.actual_capacity})

    city_link_info_dict_keys = topology_values[topology][2]

    #Store dictionaries inside oveerall dictionary to describe network topology
    for ind,key in enumerate(city_link_info_dict_keys):
        city_link_info_dict[key] = dicts_list[ind]
    
    return city_link_info_dict

def cities_link_info():

    EDFA_numbers = [1,2,3]

    #Distances can be set to 0 to observe effect of fewer links as has been done below by not using 15 links

    #Order of distances for ring - ["Lon-B","B-M","M-G","G-Le","Le-Lon"]

    ring_dists_v1 = [[81]*2,[58,55],[55] + [80]*3,[48] + [80]*3,[32] + [80]*3]
    ring_dists_v2 = [[81]*2,[58,55],[59]*5,[48] + [80]*3,[32] + [80]*3]
    ring_dists_v3 = [[40.5]*4,[58,55],[59]*5,[72]*4,[68]*4]
    ring_dists_v4 = [[40.5]*4,[58,55],[59]*5,[36]*8,[32]+[40]*6]
    #ring_dists_best_based_on_versions_tried = [[40.5]*4,[58,55],[59]*5,[36]*8,[68]*4] - very minor difference 
    #(0.003 dB) in Le-Lon edge compared to v4 (use v4 as image created for it)

    #Order of distances for mesh A - ["Lon-B","B-M","Lon-M","M-G","M-Le","G-Le","Le-Lon"]

    meshA_dists_v1 = [[40.5]*4,[58,55],[51,80,80,51],[55] + [80]*3,[59],[48] + [80]*3,[32]+[80]*3]
    meshA_dists_v2 = [[40.5]*4,[58,55],[54,52,52,52,52],[59]*5,[59],[36]*8,[32]+[40]*6]

    #Order of distances for mesh B - ["Lon-B","B-M","Lon-M","Lon-G","M-G","B-Le","M-Le","G-Le","Le-Lon"]

    meshB_dists_v1 = [[40.5]*4,[58,55],[51,80,80,51],[75] + [80]*6,[55] + [80]*3,[74]*2,[59],[48] + [80]*3,[32]+[80]*3]
    meshB_dists_v2 = [[40.5]*4,[58,55],[54,52,52,52,52],[37]*15,[59]*5,[37]*4,[59],[36]*8,[32]+[40]*6]

    #Display any effects of gain saturation on particular links and specific EDFA choices
    #Also, store dictionaries for each network topology in its respective variable
    print("Ring topology EDFA choices" + "\n" +  "\n")
    city_links_ring = network_edge_parameters(EDFA_numbers,ring_dists_v4,"ring")
    print("Mesh type A topology EDFA choices" + "\n" +  "\n")
    city_links_meshA = network_edge_parameters(EDFA_numbers,meshA_dists_v2,"meshA")
    print("Mesh type B topology EDFA choices" + "\n" +  "\n")
    city_links_meshB = network_edge_parameters(EDFA_numbers,meshB_dists_v2,"meshB")
    
    city_links_topologies = [city_links_ring,city_links_meshA,city_links_meshB]
    
    return city_links_topologies

if __name__ == "__main__":
    print("\n")
    print("*** 4B23 Coursework calculations set up ***" +  "\n")

    unformatted_output_ring,unformatted_output_meshA,unformatted_output_meshB = cities_link_info()

    #print(unformatted_output_ring)
    print("Ring topology" + "\n")
    for i in unformatted_output_ring:
        print(i + ": " + str(unformatted_output_ring[i]) + "\n")

    print("\n")
    
    #print(unformatted_output_meshA)
    print("Mesh Type A topology" + "\n")
    for i in unformatted_output_meshA:
        print(i + ": " + str(unformatted_output_meshA[i]) + "\n")

    print("\n")

    #print(unformatted_output_meshB)
    print("Mesh Type B topology" + "\n")
    for i in unformatted_output_meshB:
        print(i + ": " + str(unformatted_output_meshB[i]) + "\n")
