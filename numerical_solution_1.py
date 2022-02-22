import sys
import numpy as np
from scipy import optimize

print("""
****************************************************************
Robopool theorem for highly composite numbers:
σ_max(N)/N=∏(prime(i)/(prime(i)-1))  // ∏ - is infinite product
σ_max(N)/N ~ ln(last_prime)*eˆγ
Riemann hypothesis is true if conditions is true:
----------------------------------------------------------------
conditionos:
x>100
for π(x)=x/(ln(x)-ε), ε>1 and ε<1.1
π(x+Δ)=π(x)+1, (x+Δ)/(ln(x+Δ)-ε)-x/(ln x-ε)=1
Δ(min)=x^(1/(1+Δ(min)-1)+1)-x,
Δ>=Δ(min)
----------------------------------------------------------------
Reliability of numerical method x ϵ <100,10^14>
----------------------------------------------------------------
To end the program, press "0" and the enter.
****************************************************************
""")
# defining the necessary constants.
global epsilon
epsilon=1.00005

# enter numbers in the console.
def get_input():
    while True:
        try:
            print("Enter the number x=:")
            input_string=sys.stdin.readline()
            Num_x=int(input_string)

        except Exception:
            print("Please insert integer values")
            continue
        break
    return Num_x

# Define function delta
def func_delta(delta):
    solution=(Num_x+delta)/(np.log(Num_x+delta)-epsilon)-Num_x/(np.log(Num_x)-epsilon)-1
    return solution  

# Define function delta_min
def func_delta_min(delta_min):
    solution=Num_x**(1/(Num_x+delta_min-1)+1)-Num_x-delta_min
    return solution

# Infinite while loop console. Main program
while True:
    Num_x=get_input()   
    # end of program    
    if Num_x == 0:
        break         
    # initial function delta
    initial_num=0
    delta_solution=optimize.broyden1(func_delta,initial_num )
    # initial function delta_min
    delta_min_solution=optimize.broyden1(func_delta_min,initial_num)
    # listing results
    residual=(Num_x+delta_solution)/(np.log(Num_x+delta_solution)-epsilon)-Num_x/(np.log(Num_x)-epsilon)-1
    print(f'the residual of the numerical calculation Δ: {residual}')
    print(f'Δ is: {delta_solution}')
    print(f'Δ(min) is: {delta_min_solution}')
    if delta_solution>=delta_min_solution:
        print("delta_solution>=delta_min_solution")
    else:
        print("the condition failed")
       
    
