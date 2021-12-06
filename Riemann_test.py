# import library
import sys
import math

print("************************************************************************************************************************") 
print("Riemann hypothesis")
print("")
print("Author: Ing. Robert Polak")
print("website: https://www.robopol.sk")
print("Purpose:")
print("In 1984 Guy Robin showed that the Riemann Hypothesis is true if and only if:")
print("σ(n)/(n*ln(ln(n))) < eˆγ")
print("eˆγ = 1.781072417990197985236504103107...")
print("The Riemann Hypothesis Says 5040 is the Last:")
print("σ(5040)/(5040*ln(ln(5040))) = 1.7909733665348811333619013505910…")           
print("")
print("Ideal pattern sequence that maximizes the Guy Robin index")
print("Guy Robin calculation procedure for sequence: N=p(1)ˆa*p(2)ˆb*p(3)ˆc...p(i)ˆ1; p(i)=2,3,5,7 ...; p(i) ϵ prime numbers,")
print("")
print("Copyright notice: This code is Open Source and should remain so.")
print("")
print("To end the program, press 0 and the enter.")           
print("************************************************************************************************************************")

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

#  function ideal numbers
def basic_sequence(n):
    # defining the necessary constants   
    k=1; sigma_seq_1=12; numb_seq_1=6; index=2
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
            field_primes_count_n.append(divis_1)
            numb_seq_1*=divis_1
            sigma_seq_1*=(divis_1+1)
            index+=1
        if vysledok2==1:
            field_primes_count_n.append(divis_2)
            numb_seq_1*=divis_2
            sigma_seq_1*=(divis_2+1)
            index+=1
        k+=1
    # calculate Guy Robin index
    if sigma_seq_1 >=10**300: 
        logarithm=math.log(math.log(numb_seq_1))       
        divisor=multipl(logarithm,numb_seq_1)
        guy_robin=longDivision(sigma_seq_1,divisor)                
    else:
        guy_robin=sigma_seq_1/(numb_seq_1*math.log(math.log(numb_seq_1)))        
    
    # computing for ideal number
    constant_1=0;constant_2=0;guy_max=guy_robin;por_cislo=0;cislo=2;index_1=0
    sigma_nova=sigma_seq_1; sigma_pomocna=1
    numb_seq_3=numb_seq_1
    while constant_2<2:
        sigma_pomocna=1
        for field in field_primes_count_n:
            if field!=cislo: sigma_pomocna*=(field+1)                        
        mocnina=1
        pomocne_cislo=field_primes_count_n[index_1]
        while constant_1<1:
            mocnina+=1
            sigma_nova+=(cislo**mocnina)*sigma_pomocna
            numb_seq_3*=cislo
            pomocne_cislo+=cislo**mocnina
            # calculate Guy Robin index
            if sigma_nova >=10**300: 
                logarithm=math.log(math.log(numb_seq_3))       
                divisor=multipl(logarithm,numb_seq_3)
                guy_robin=longDivision(sigma_nova,divisor)                
            else:
                guy_robin=sigma_nova/(numb_seq_3*math.log(math.log(numb_seq_3)))
            # conditions for guy_max
            if guy_robin>guy_max:
                constant_1=0; constant_2=0
                guy_max=guy_robin
                num_max=numb_seq_3
            else:
                constant_1=1; constant_2+=1
        # field change
        del field_primes_count_n[index_1]
        field_primes_count_n.insert(index_1,pomocne_cislo)    
        por_cislo+=1
        index_1+=1        
        cislo=field_primes_count_n[index_1]                
        constant_1=0        
    # conditions for the statement
    if num_max>10**10000:
        print(f'The number has: {len(str(num_max))} digits')
        print(f'Guy Robin is: {guy_max}')
    else:     
        print(f'Number is: {num_max}')
        print(f'Guy Robin is: {guy_max}')    
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
    # initialization function ideal number
    basic_sequence(numb)       