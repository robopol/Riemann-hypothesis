# ************************************************************************************************************************************
# Riemann hypothesis
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
#
#           1) Guy Robin calculation procedure for sequence: N=p(1)*p(2)*p(3)...p(i); p(i)=2,3,5,7 ...; p(i) ϵ prime numbers
#           2) Guy Robin calculation procedure for sequence: N=p(1)ˆa*p(2)ˆb*p(3)ˆc...p(i)ˆ1; p(i)=2,3,5,7 ...; p(i) ϵ prime numbers,
#              a,b,c ϵ N; a>b>c; example N=2*3*4*5*6...
#           3) Ideal pattern sequence based on 1) and 2) that maximizes the Guy robin index
# 'type: console program' 
# Copyright notice: This code is Open Source and should remain so.
# ************************************************************************************************************************************
import sys
import math
from itertools import combinations
import matplotlib.pyplot as plt

print("**************************************************************************")
print("Riemann hypothesis:")
print("")
print("To select a sequence, press 0,1,2, 3 -(it is not available yet)")
print("0 - to end the program, press 0 and the enter.")
print("1 - infinite prime number sequence: N=2*3*5*7*11*13...")
print("2 - infinite sequence: N=2*3*4*5*6*7*8...")
print("3 - ideal pattern: it is not available yet")
print("To end the program, press 0 and the enter.")
print("**************************************************************************")

# defining the necessary constants.
e_gama=1.7909733665348811333619013505910

# enter numbers in the console.
def get_input():
    while True:
        try:
            input_string=sys.stdin.readline()
            number_int=int(input_string)
        except Exception:
            print("Please insert only integer values")
            continue
        break
    return number_int

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
        # computing decemposition        
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
    sigma_value=1; w=1    
    for sig_i in sigma_field:
        auxiliary_field=list(combinations(sigma_field,w))        
        delta_sigma=1
        for aux_field_i in auxiliary_field:
            delta_sigma=1
            for x in aux_field_i:
                delta_sigma=x*delta_sigma
            sigma_value+=delta_sigma        
        w+=1   
    return sigma_value

#  function sequence (1)
def sequence_1(count_n):
    # creation of a field of primes
    k=1; divis_2=1
    field_primes_count_n=[2,3]
    while len(field_primes_count_n) <= count_n:
        divis_1=6*k-1
        divis_2=6*k+1
        k+=1; t=0; r=0
        for i in field_primes_count_n:                          
            if divis_1 % i !=0:
                t+=1
                if t==len(field_primes_count_n):
                    field_primes_count_n.append(divis_1)
            if divis_2 % i !=0:
                r+=1
                if r==len(field_primes_count_n):
                    field_primes_count_n.append(divis_2)
    
    field_sequence_1=[]; number=1
    for j in field_primes_count_n:
        number*=j
        if number >2:
            # call decomposition function
            decomposition=decomp_number(number)
            # call sigma function
            sigma_value=sigma(decomposition)
            # calculate Guy Robin index for sequence
            guy_robin=sigma_value/(number*math.log(math.log(number)))
            # field sequence 1
            field_sequence_1.append([number,guy_robin])    
    return field_sequence_1

#  function sequence (2)
def sequence_2(count_n):
    # creation of a field sequence 2
    field_sequence_2=[]; number=1
    for i in range(1,count_n+3):
        number*=i
        if number >2:
            # call decomposition function
            decomposition=decomp_number(number)
            # call sigma function
            sigma_value=sigma(decomposition)
            # calculate Guy Robin index for sequence
            guy_robin=sigma_value/(number*math.log(math.log(number)))
            # field sequence 2
            field_sequence_2.append([number,guy_robin])
    return field_sequence_2

# Infinite while loop console. Main program
while True:
    print("Enter the mod (0,1,2):")
    sequen=get_input()   
    # end of program (0)   
    if sequen == 0:
        break     
    # sequence (1)
    if sequen==1:
        print("Enter the number of members in sequence '1', max 20 (memory overflow):")
        count_n=get_input()
        field_sequence_1=sequence_1(count_n)
        print("--wait for the result--")
        print("Field sequence '1' is:")
        print(field_sequence_1)
        index_i=1; field_chart_1=[]
        for field_1_i in field_sequence_1:
            field_chart_1.append(field_1_i[-1])           
        # plot a sequence graph (1)
        plt.title(" Chart Guy Robin sequence no.'1'")
        plt.grid(True)
        plt.xlabel("serial number in sequence of numbers N =6,30,210 ...")
        plt.ylabel("value Guy - Robin index")
        plt.plot(field_chart_1) 
        plt.show()    

    # sequence (2)
    if sequen==2:
        print("Enter the number of members in sequence '2', max 100 (memory overflow):")
        count_n=get_input()
        field_sequence_2=sequence_2(count_n)
        print("--wait for the result--")
        print("Field sequence '2' is:")
        print(field_sequence_2)
        index_i=1; field_chart_2=[]
        for field_1_i in field_sequence_2:
            field_chart_2.append(field_1_i[-1])            
        # plot a sequence graph (2)
        plt.title(" Chart Guy Robin sequence no.'2'")
        plt.grid(True)
        plt.xlabel("serial number in sequence of numbers N =6,24,120 ...")
        plt.ylabel("value Guy - Robin index")
        plt.plot(field_chart_2) 
        plt.show()