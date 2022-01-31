# **************************************************************************************
# σ(max) for highly composite numbers
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
# 'type: console program' 
# Copyright notice: This code is Open Source and should remain so.
# **************************************************************************************
import sys
import math

print("*****************************************************************")
print("Robopool theorem for highly composite numbers:")
print("σ_max(N)/N=∏(prime(i)/(prime(i)-1))  // ∏ - is infinite product")
print("σ_max(N)/N ~ ln(last_prime)*eˆγ")
print("")
print("To end the program, press 0 and the enter.")
print("*****************************************************************")

# defining the necessary constants.
e_gama=1.7909733665348811333619013505910
global field_primes_count_n
global count_n

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

#  function basic sequence primes
def basic_sequence(n):
    # defining the necessary constants   
    k=1; sigma_m=2; cislo=6; index=2
    field_primes_count_n=[2,3]
    while index < n:
        divis_1=6*k-1
        divis_2=6*k+1
        vysledok1=1;vysledok2=1
        hranica=round(math.sqrt(divis_2)+1)        
        for field in field_primes_count_n:
            if divis_1 % field==0:vysledok1=0
            if divis_2 % field==0:vysledok2=0    
            if field > hranica:
                break                
        if vysledok1==1:
            last_prime=divis_1
            field_primes_count_n.append(divis_1)
            cislo*=divis_1
            sigma_m*=(divis_1-1)
            index+=1
        if vysledok2==1:
            last_prime=divis_2
            field_primes_count_n.append(divis_2)
            cislo*=divis_2
            sigma_m*=(divis_2-1)
            index+=1
        k+=1                    
    # conditions for the statement
    # last prime
    print(f'The last prime is: {last_prime}')
    # number
    if cislo>10**10000:
        print(f'The number has: {len(str(cislo))} digits')        
    else:     
        print(f'Number is: {cislo}')
             
    # calculate Guy Robin index
    if cislo >=10**300:
        log_lastprime=math.log(last_prime)      
        ratio=longDivision(cislo,sigma_m)               
        print(f'log(last_prime) is: {log_lastprime}')        
        print(f'σ_max(n)/N is: {ratio}')
        print(f'(eˆγ*log(last_prime) is: {e_gama*log_lastprime}')        
                
    else:
        log_lastprime=math.log(last_prime)
        print(f'log(last_prime) is: {log_lastprime}')        
        print(f'σ_max(n)/N is: {cislo/sigma_m}')
        print(f'(eˆγ*log(last_prime) is: {e_gama*log_lastprime}')        
        
    return     

# function divide large numbers
def longDivision(number, divisor):
    ans=str(number//divisor)
    remainder=1
    if number % divisor !=0:
        ans=ans+"."
        remainder=number % divisor    
    for i in range(1,11):
        temp=remainder*10
        ans+=str(temp//divisor)
        remainder=temp % divisor
    ans=float(ans)
    return ans

# function: multiplication large numbers
def multipl(float_num, number):
    summ_temp=0; point=0; int_numb=""
    divis=0
    temp=str(float_num)
    for i in range(0,len(temp)-1):        
        if temp[i]==".":
            point=i
            break
        else:
            int_numb+=temp[i]
    for j in range(point+1,len(temp)):
        k=(int(temp[j])*int(number))//10**(j-point)
        summ_temp+=int(k)    
    if int(temp[0])==0: divis=summ_temp                
    if int(temp[0])>0: divis=number*int(int_numb)+summ_temp        
    return divis

# Infinite while loop console. Main program
while True:
    numb=get_input()   
    # end of program    
    if numb == 0:
        break         
    # initialization of the number decomposition function
    basic_sequence(numb)        