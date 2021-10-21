# import library
import math
import matplotlib.pyplot as plt
from itertools import combinations
from tkinter import *
from tkinter import ttk
from tkinter import messagebox
from numpy import integer

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
print("1) Guy Robin calculation procedure for sequence: N=p(1)*p(2)*p(3)...p(i); p(i)=2,3,5,7 ...; p(i) ϵ prime numbers")
print("2) Guy Robin calculation procedure for sequence: N=p(1)ˆa*p(2)ˆb*p(3)ˆc...p(i)ˆ1; p(i)=2,3,5,7 ...; p(i) ϵ prime numbers,")           
print("   a,b,c ϵ N; a>b>c; example N=2*3*4*5*6...")
print("3) Ideal pattern sequence based on 1) and 2) that maximizes the Guy Robin index")
print("")
print("Copyright notice: This code is Open Source and should remain so.")           
print("************************************************************************************************************************") 

# Main program - tkinker
root=Tk()
root.title("Riemann hypothesis, version:2.0")

# tkinker forms  
frame_riemann=LabelFrame(root, text="Calculate Guy Robin", padx=20, pady=20)
frame_riemann.pack(padx=20, pady=20)
e=Entry(frame_riemann,width=20,borderwidth=2)
e.pack()
r=IntVar()

# defining the necessary constants.
e_gama=1.7909733665348811333619013505910

def get_input():    
    try:
        input_string=e.get()
        number_int=int(input_string)
    except Exception:            
        messagebox.showwarning(title=Warning, message="Please insert only integer values")       
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

# function calculate sigma_iter
def sigma_comb(iterable, r):    
    pool = tuple(iterable)
    delta_sigma=1
    for k in range(0,r):
        delta_sigma *=pool[k]    
    n = len(pool)    
    indices = list(range(r))        
    while True:        
        a=1
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return delta_sigma
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1        
        for i in indices:
            a*=int(pool[i])         
        delta_sigma +=a           

# function calculate sigma
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
    if len(sigma_field)>27:
        for sig_i in sigma_field:
            delta_sigma=sigma_comb(sigma_field,w)        
            sigma_value+=delta_sigma
            w+=1           
        return sigma_value
    else:        
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

# function sequence (3)
def sequence_3(count_n):
    # creation of a field sequence 3
    field_sequence_3=[] 
    # creation of a field primes
    field_primes_init=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307]
    # creating a temporary field, defining the necessary constants
    field_temp_min=[2,1,1,1] 
    guy_robin_actual=1; index=1
    num_cycle=0
    # main cycle for count_n
    while num_cycle < count_n:        
        index=1        
        while index <= 2:
            if index==1:
                # cyklus 1
                i_1=0; i_2=0; field_temp=list(field_temp_min); guy_robin_temp=0.2            
                while i_2 <2:            
                    k=0; number=1; decomposition=[]
                    for field_i in field_temp:
                        number*=field_primes_init[k]**field_i
                        k+=1                    
                    # call decomposition function    
                    decomposition=decomp_number(number)
                    # call sigma function
                    sigma_value=sigma(decomposition)
                    # calculate Guy Robin index for sequence
                    guy_robin_actual=sigma_value/(number*math.log(math.log(number)))                
                    if guy_robin_actual >=1.72:
                        y=0
                        for field_3_i in field_sequence_3:
                            if field_3_i == [number,guy_robin_actual]: y+=1
                        if y==0:
                            field_sequence_3.append([number,guy_robin_actual])
                            num_cycle+=1                
                    # condition for the next index
                    if guy_robin_actual > guy_robin_temp:
                        guy_robin_temp=guy_robin_actual
                        field_temp[i_1]+=1; i_2=0
                        if i_1==0: a=field_temp[0]
                    else:
                        if i_2==0: field_temp[i_1]-=1           
                        if i_1<len(field_temp):
                            if i_2==1: field_temp[i_1]-=1
                            i_1+=1
                            field_temp[i_1]+=1
                            i_2+=1
                            if i_2==2: field_temp[i_1]-=1
                        else: break                                    
                # next step
                index=2
                # defining new field_temp_min
                field_next=[]
                for i in range(0,len(field_temp)):
                    if field_temp[i]-1<=1:
                        field_next.append(1)
                    else:
                        field_next.append(field_temp[i]-1)
                field_next.append(1)                            
            if index==2:
                # cyklus 2
                i_1=0; i_2=0; field_temp=list(field_temp_min); guy_robin_temp=0.2
                field_temp[0]=a-1            
                while i_2 <2:            
                    k=0; number=1; decomposition=[]
                    for field_i in field_temp:
                        number*=field_primes_init[k]**field_i
                        k+=1                    
                    # call decomposition function    
                    decomposition=decomp_number(number)
                    # call sigma function
                    sigma_value=sigma(decomposition)
                    # calculate Guy Robin index for sequence
                    guy_robin_actual=sigma_value/(number*math.log(math.log(number)))
                    if guy_robin_actual >=1.72:
                        y=0
                        for field_3_i in field_sequence_3:
                            if field_3_i == [number,guy_robin_actual]: y+=1
                        if y==0:
                            field_sequence_3.append([number,guy_robin_actual])
                            num_cycle+=1
                    # condition for the next index
                    if guy_robin_actual > guy_robin_temp:
                        guy_robin_temp=guy_robin_actual
                        field_temp[i_1]+=1; i_2=0
                    else:
                        if i_1 >0:
                            if i_2==0:field_temp[i_1]-=1           
                            if i_1<len(field_temp):
                                if i_2==1: field_temp[i_1]-=1
                                i_1+=1
                                field_temp[i_1]+=1
                                i_2+=1
                                if i_2==2: field_temp[i_1]-=1
                            else: break                                                        
                        else:
                            i_1=1; field_temp[i_1]+=1; guy_robin_temp=0.2                        
                # next step
                index=3                
                # field_temp_min initialization
                field_temp_min=list(field_next)              
    # field_sequence_3 arrangement
    field_sequence_3.sort()
    return field_sequence_3

# function sequence()
def sequence():    
    answer=True
    # for sequence (1)    
    if r.get()==1:        
        if e.get()=="":
            messagebox.showwarning(title=Warning, message="Enter the number of members in sequence '1'-(only integer value):")
        count_n=get_input()        
        print("--wait for the result--")
        field_sequence_1=sequence_1(count_n)        
        print("Field sequence '1' is:")
        print(field_sequence_1)
        index_i=1; field_chart_1=[]
        for field_1_i in field_sequence_1:
            field_chart_1.append(field_1_i[-1])           
        # plot a sequence graph (1)
        plt.title("Chart Guy Robin sequence no.'1'")
        plt.grid(True)
        plt.xlabel("x number in sequence of numbers N =6,30,210 ...")
        plt.ylabel("value Guy Robin index")
        plt.plot(field_chart_1) 
        plt.show()
    # for sequence (2)
    if r.get()==2:        
        if e.get()=="":
            messagebox.showwarning(title=Warning, message="Enter the number of members in sequence '2'-(only integer value):")
        count_n=get_input()        
        print("--wait for the result--")
        field_sequence_2=sequence_2(count_n)        
        print("Field sequence '2' is:")
        print(field_sequence_2)
        index_i=1; field_chart_2=[]
        for field_1_i in field_sequence_2:
            field_chart_2.append(field_1_i[-1])            
        # plot a sequence graph (2)
        plt.title(" Chart Guy Robin sequence no.'2'")
        plt.grid(True)
        plt.xlabel("x number in sequence of numbers N =6,24,120 ...")
        plt.ylabel("value Guy Robin index")
        plt.plot(field_chart_2) 
        plt.show()
    # for sequence (3)
    if r.get()==3:        
        if e.get()=="":
            messagebox.showwarning(title=Warning, message="Enter the number of members in sequence '3'-(only integer value):")
        count_n=get_input()
        print("--wait for the result--")
        field_sequence_3=sequence_3(count_n)
        print("Field sequence '3' is:")
        print(field_sequence_3)
        index_i=1; field_chart_3=[]
        for field_3_i in field_sequence_3:
            field_chart_3.append(field_3_i[-1])
        # plot a sequence graph (3)
        plt.title(" Chart Guy Robin sequence no.'3'- ideal numbers")
        plt.grid(True)
        plt.xlabel("x number in sequence of ideal numbers")
        plt.ylabel("value Guy Robin index")
        plt.plot(field_chart_3) 
        plt.show()
    return

# tkinter forms
b_start=Button(frame_riemann, text="Start",padx=10, pady=5, command=sequence)
myLabel=Label(frame_riemann,text="Enter the number of members in sequence",)
myLabel.pack(padx=10, pady=5)

Radiobutton(frame_riemann,text="1 - infinite prime number sequence: N=2*3*5*7*11*13...",variable=r, value=1).pack(anchor=W)
Radiobutton(frame_riemann,text="2 - infinite sequence: N=2*3*4*5*6*7*8...",variable=r, value=2).pack(anchor=W)
Radiobutton(frame_riemann,text="3 - ideal pattern sequence: that maximizes the Guy Robin index",variable=r, value=3).pack(anchor=W)
b_start.pack(padx=10, pady=5)


root.mainloop()

