# **************************************************************************************
# Guy Robin index
#
# Author: ing. Robert Polak
# Contact Info: robopol@robopol.sk
# website: https://www.robopol.sk
# Purpose:
#           In 1984 Guy Robin showed that the Riemann Hypothesis is true if and only if:
#           σ(n)/(n*ln(ln(n))) < eˆγ
#           eˆγ = 1.781072417990197985236504103107...
#           The Riemann Hypothesis Says 5040 is the Last:
#           σ(5040)/(5040*ln(ln(5040))) = 1.7909733665348811333619013505910…
# 'type: console program' 
# Copyright notice: This code is Open Source and should remain so.
# **************************************************************************************
import sys
import math
# import library for all combinations, to get started, install the library: pip install more-itertools
from more_itertools import powerset

print("*****************************************************************")
print("Guy Robin index calculation:")
print("")
print("To end the program, press 0 and the enter.")
print("*****************************************************************")

# defining the necessary constants.
e_gama=1.7909733665348811333619013505910

# enter numbers in the console.
def get_input():
    while True:
        try:
            print("Enter the number:")
            input_string=sys.stdin.readline()
            number_int=int(input_string)
        except Exception:
            print("Please insert only integer values")
            continue
        break
    return number_int

# Filtering prime, method creator: robopol.sk
def prime_fun(numb):
    # defining the necessary constants.
    b_first=251; a=2
    cycle_end=cycle= int(a**b_first % numb)
    complet_end=complet=int(b_first)
    euler=(numb - 1)//2
    k=euler
    # We define an array of residues, numbers.
    field_num=[b_first]
    field_res=[cycle]
    # We will fill these fields.
    while complet < euler:
        if 2*complet >= euler:
            break
        cycle=cycle**2
        if cycle > numb:
            cycle=cycle % numb
            if cycle == 0: cycle=numb
        complet=2*complet
        field_num=field_num + [complet]
        field_res=field_res + [cycle]
    # Main algorithm.                       
    while k > b_first:
        index=0
        for num in field_num:
            if num < k: 
                maximum = num
                max_cycle = field_res [index]
            index +=1
        complet_end += maximum
        cycle_end = cycle_end*max_cycle
        if cycle_end > numb:
            cycle_end=cycle_end % numb
            if cycle_end == 0: cycle_end=numb
        k= euler - complet_end
    remainder = int((cycle_end * a ** k) % numb)
    if numb-remainder == 1:
        remainder=-1
    return remainder 

# number decomposition function
def decomp_number(n):
    # defining the necessary constants.
    big_num=int(10**13)
    prime_field=[2,3]
    decomp_field_basic=[]
    divisor_1=divisor_2=1
    k=0; rest=int(n)
    odpocet_1=0; odpocet_2=0 
    # decomposition function for 2,3
    for prime in prime_field:
        w=0
        while rest % prime == 0:
            w+=1
            rest=rest//prime            
        if w>=1:
            decomp_field_basic.append([prime,w])                            
    # decomposition function for >3...
    range_big=round(math.sqrt(big_num)+1)
    while rest > 1:     
        k+=1;i=0                
        divisor_1=6*k-1
        divisor_2=6*k+1
        # computing for rest<big_num
        if rest<=big_num:
            if divisor_1 <= range_big:
                while rest % divisor_1 == 0:
                    i+=1
                    rest=rest//divisor_1
                if i >=1:
                    decomp_field_basic.append([divisor_1,i])
                i=0
                while rest % divisor_2 == 0:
                    i+=1
                    rest=rest//divisor_2
                if i >=1:
                    decomp_field_basic.append([divisor_2,i])
            else:
                decomp_field_basic.append([rest,1])
                break
        # computing for rest>big_num
        if rest>big_num:
            if odpocet_1 == 0:
                # Call functions to find a prime number
                index= prime_fun(rest)                       
                if index ==1 or index ==-1:
                    decomp_field_basic.append([rest,1])
                    break
                else:                                           
                    while rest % divisor_1 == 0:
                        i+=1
                        rest=rest//divisor_1
                        odpocet_1=0
                    if i >=1:
                        decomp_field_basic.append([divisor_1,i])
                    i=0
                    while rest % divisor_2 == 0:
                        i+=1
                        rest=rest//divisor_2
                        odpocet_1=0
                    if i >=1:
                            decomp_field_basic.append([divisor_2,i])                    
            else:
                while rest % divisor_1 == 0:
                    i+=1
                    rest=rest//divisor_1
                    odpocet_1=0
                if i >=1:
                    decomp_field_basic.append([divisor_1,i])
                i=0
                while rest % divisor_2 == 0:
                    i+=1
                    rest=rest//divisor_2
                    odpocet_1=0
                if i >=1:
                    decomp_field_basic.append([divisor_2,i])
    return decomp_field_basic

# Function calculate sigma
def sigma(decomposition):
    # defining the necessary constants.
    sigma_field=[]    
    for decom in decomposition:
        sigma_decom_value=0        
        for i in range(1,decom[1]+1):
            sigma_decom_value +=int(decom[0]**(i))
        sigma_field.append(sigma_decom_value)
    # calculate sigma
    sigma_value=0   
    sigma_power=list(powerset(sigma_field))
    sigma_power[0]=(1,)    
    for sigma_power_i in sigma_power:
        delta_sigma=1; j=0        
        for x in sigma_power_i:
            delta_sigma=x*delta_sigma                    
        sigma_value+=delta_sigma    
    print(f'sigma is: {sigma_value}')
    return sigma_value

# Infinite while loop console. Main program
while True:
    number=get_input()   
    # end of program    
    if number == 0:
        break 
    # initialization of the number decomposition function
    decomposition=decomp_number(number)      
    print("Number decomposition is:")    
    v=1
    for numb in decomposition:
        u=1                        
        for i in numb:
            if u == len(numb):
                print (f'{i}',end='')
            else:
                u+=1
                print (f'{i}ˆ',end='')                              
        if v < len(decomposition):
            print(f' * ',end='')
            v+=1      
    print("")
    # initialization function calculate sigma
    print("--wait for the result--")
    sigma_value=sigma(decomposition)
    # calculate Guy Robin index
    guy_robin=sigma_value/(number*math.log(math.log(number)))
    print(f'Guy Robin index is: {guy_robin}')       