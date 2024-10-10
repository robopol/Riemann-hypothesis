# Import libraries
import math
import matplotlib.pyplot as plt
from tkinter import *
from tkinter import ttk
from tkinter import messagebox
import threading
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
matplotlib.use('TkAgg')

# Display introductory text in console (optional, can be removed if not needed)
print("************************************************************************************************************************") 
print("Riemann hypothesis")
print("")
print("Author: Ing. Robert Polak")
print("Website: https://www.robopol.sk")
print("Purpose:")
print("In 1984 Guy Robin showed that the Riemann Hypothesis is true if and only if:")
print("σ(n)/(n*ln(ln(n))) < e^γ")
print("e^γ = 1.781072417990197985236504103107...")
print("The Riemann Hypothesis Says 5040 is the Last:")
print("σ(5040)/(5040*ln(ln(5040))) = 1.7909733665348811333619013505910…")           
print("1) Guy Robin calculation procedure for sequence: N=p(1)*p(2)*p(3)...p(i); p(i)=2,3,5,7 ...; p(i) ∈ prime numbers")
print("2) Guy Robin calculation procedure for sequence: N=p(1)^a*p(2)^b*p(3)^c...p(i)^1; p(i)=2,3,5,7 ...; p(i) ∈ prime numbers,")           
print("   a,b,c ∈ N; a>b>c; example N=2*3*4*5*6...")
print("3) Ideal pattern sequence based on 1) and 2) that maximizes the Guy Robin index")
print("")
print("Copyright notice: This code is Open Source and should remain so.")           
print("************************************************************************************************************************") 

# Main program - Tkinter
root = Tk()
root.title("Riemann Hypothesis, version:2.0")
root.geometry("900x700")
root.configure(bg='#2b2b2b')  # Set dark background color

# Apply a dark theme using ttk styles
style = ttk.Style()
style.theme_use('clam')
style.configure('.', background='#2b2b2b', foreground='white', fieldbackground='#2b2b2b')
style.configure('TLabel', background='#2b2b2b', foreground='white')
style.configure('TButton', background='#444444', foreground='white')
style.configure('TEntry', fieldbackground='#444444', foreground='white')
style.configure('TRadiobutton', background='#2b2b2b', foreground='white')
style.configure('TFrame', background='#2b2b2b')

# Tkinter frames
frame_riemann = ttk.LabelFrame(root, text="Calculate Guy Robin", padding=(20, 20))
frame_riemann.pack(padx=20, pady=10, fill=X)

input_frame = ttk.Frame(frame_riemann)
input_frame.pack(pady=10, fill=X)

# Entry for number of members
entry_label = ttk.Label(input_frame, text="Enter the number of members in sequence:")
entry_label.pack(anchor=W, pady=(0, 5))
e = ttk.Entry(input_frame, width=20)
e.pack(anchor=W, pady=(0, 10))

# Radio buttons for sequence selection
r = IntVar()
sequence_frame = ttk.LabelFrame(frame_riemann, text="Select Sequence", padding=10)
sequence_frame.pack(fill=X, pady=10)

radiobutton1 = ttk.Radiobutton(sequence_frame, text="1 - Infinite prime number sequence: N=2*3*5*7*11*13...", variable=r, value=1)
radiobutton1.pack(anchor=W, pady=2)
radiobutton2 = ttk.Radiobutton(sequence_frame, text="2 - Infinite sequence: N=2*3*4*5*6*7*8...", variable=r, value=2)
radiobutton2.pack(anchor=W, pady=2)
radiobutton3 = ttk.Radiobutton(sequence_frame, text="3 - Ideal pattern sequence: that maximizes the Guy Robin index", variable=r, value=3)
radiobutton3.pack(anchor=W, pady=2)

# Start button
start_button = ttk.Button(frame_riemann, text="Start", command=lambda: run_sequence())
start_button.pack(pady=10)

# Progress label
progress_label = ttk.Label(frame_riemann, text="", background='#2b2b2b', foreground='yellow')
progress_label.pack()

# Output frame
output_frame = ttk.Frame(root)
output_frame.pack(padx=20, pady=10, fill=BOTH, expand=True)

# Defining the necessary constants.
e_gama = 1.7909733665348811333619013505910
global field_primes_count_n

# Function to get user input
def get_input():    
    try:
        input_string = e.get()
        number_int = int(input_string)
    except Exception:            
        messagebox.showwarning(title="Warning", message="Please insert only integer values")
        number_int = None
    return number_int

# Number decomposition function
def decomp_number(n):
    big_num = int(10**13)
    prime_field = [2, 3]
    decomp_field_basic = []
    divisor_1 = divisor_2 = 1
    k = 0
    rest = int(n)
    odpocet_1 = 0
    odpocet_2 = 0 
    # Decomposition for 2,3
    for prime in prime_field:
        w = 0
        while rest % prime == 0:
            w += 1
            rest = rest // prime            
        if w >= 1:
            decomp_field_basic.append([prime, w])                            
    # Decomposition for >3...
    range_big = round(math.sqrt(big_num) + 1)
    while rest > 1:     
        k += 1
        i = 0                
        divisor_1 = 6 * k - 1
        divisor_2 = 6 * k + 1
        # Computing decomposition        
        if divisor_1 <= range_big:
            while rest % divisor_1 == 0:
                i += 1
                rest = rest // divisor_1
            if i >=1:
                decomp_field_basic.append([divisor_1, i])
            i = 0
            while rest % divisor_2 == 0:
                i += 1
                rest = rest // divisor_2
            if i >=1:
                decomp_field_basic.append([divisor_2, i])
        else:
            decomp_field_basic.append([rest,1])
            break        
    return decomp_field_basic         

# Function to divide large numbers
def longDivision(numb, divisor):
    ans = str(numb // divisor)
    remainder = 1    
    if numb % divisor != 0:
        ans = ans + "."
        remainder = numb % divisor    
    for i in range(1, 11):
        temp = remainder * 10
        ans += str(temp // divisor)
        remainder = temp % divisor
    ans = float(ans)
    return ans

# Function: multiplication large numbers
def multipl(float_num, number):
    point = 0
    int_numb = ""
    summ_temp = 0
    multiple_sum = 0
    temp = str(float_num)
    for i in range(0, len(temp) - 1):        
        if temp[i] == ".":
            point = i
            break
        else:
            int_numb += temp[i]    
    for j in range(point + 1, len(temp)):
        k = (int(temp[j]) * int(number)) // 10 ** (j - point)
        summ_temp += int(k)            
    if int(temp[0]) == 0:
        multiple_sum = summ_temp                
    if int(temp[0]) > 0:
        multiple_sum = number * int(int_numb) + summ_temp            
    return multiple_sum

# New function to calculate sigma, Robopol's theorem
def sigma_new(decomposition):
    sigma_last = 1
    decomp_correction = []     
    for decom in decomposition:
        number_correction = 0
        for i in range(1, decom[1] + 1):
            number_correction += int(decom[0] ** (i))
        decomp_correction.append(number_correction)
    # Calculate sigma
    for numb in decomp_correction:
        sigma_last *= int(numb + 1)
    return sigma_last       

# Function to create sequence of primes
def basic_sequence(count_n):
    k = 1
    index = 2
    field_primes_count_n = [2, 3]    
    while index < count_n:
        divis_1 = 6 * k - 1
        divis_2 = 6 * k + 1
        const_1 = 1
        const_2 = 1
        border = round(math.sqrt(divis_2) + 1)
        for i in field_primes_count_n:
            if divis_1 % i == 0:
                const_1 = 0
            if divis_2 % i == 0:
                const_2 = 0    
            if i > border:
                break                          
        if const_1 == 1: 
            field_primes_count_n.append(divis_1)
            index += 1
        if const_2 == 1: 
            field_primes_count_n.append(divis_2)
            index += 1   
        k += 1                        
    return field_primes_count_n

# Function for sequence (1)
def sequence_1(count_n):
    field_sequence_1 = []
    field_primes_count_n = basic_sequence(count_n)
    number = 1
    sigma_value = 12
    for j in field_primes_count_n:
        number *= j
        if number > 6:            
            # New sigma function, Robopol's theorem
            sigma_value *= (j + 1)
            # Calculate Guy Robin index for sequence
            if sigma_value > 10 ** 300:
                logarithm = math.log(math.log(number))       
                divis = multipl(logarithm, number)                
                guy_robin = longDivision(sigma_value, divis)
                # Field sequence 1
                field_sequence_1.append([number, guy_robin])
            else:
                guy_robin = sigma_value / (number * math.log(math.log(number)))            
                # Field sequence 1        
                field_sequence_1.append([number, guy_robin])                  
    return field_sequence_1

# Function for sequence (2)
def sequence_2(count_n):
    field_sequence_2 = []
    number = 1
    sigma_value = 3
    for i in range(1, count_n + 3):
        number *= i
        if number > 2:
            # Call decomposition function
            decomposition = decomp_number(number)
            # Call new sigma function, Robopol's theorem
            sigma_value = sigma_new(decomposition)
            # Calculate Guy Robin index for sequence
            if sigma_value >= 10 ** 300:
                logarithm = math.log(math.log(number))       
                divis = multipl(logarithm, number)
                guy_robin = longDivision(sigma_value, divis)
                # Field sequence 2
                field_sequence_2.append([number, guy_robin])
            else:
                guy_robin = sigma_value / (number * math.log(math.log(number)))
                # Field sequence 2
                field_sequence_2.append([number, guy_robin])            
    return field_sequence_2

# Function for sequence (3)
def sequence_3(count_n):
    field_sequence_3 = [] 
    field_primes_init = basic_sequence(count_n)
    field_temp_min = [2, 1, 1, 1] 
    guy_robin_actual = 1
    index = 1
    num_cycle = 0
    while num_cycle < count_n:        
        index = 1        
        while index <= 2:
            if index == 1:
                # Cycle 1
                i_1 = 0
                i_2 = 0
                field_temp = list(field_temp_min)
                guy_robin_temp = 0.2            
                while i_2 < 2 and num_cycle < count_n:            
                    k = 0
                    number = 1
                    decomposition = []
                    for field_i in field_temp:
                        number *= field_primes_init[k] ** field_i
                        k += 1                    
                    # Call decomposition function    
                    decomposition = decomp_number(number)
                    # Call new sigma function, Robopol's theorem
                    sigma_value = sigma_new(decomposition)
                    # Calculate Guy Robin index for sequence
                    if sigma_value > 10 ** 300:
                        logarithm = math.log(math.log(number))       
                        divis = multipl(logarithm, number)                
                        guy_robin_actual = longDivision(sigma_value, divis)                
                    else:
                        guy_robin_actual = sigma_value / (number * math.log(math.log(number)))              
                    if guy_robin_actual >= 1.72:
                        y = 0
                        for field_3_i in field_sequence_3:
                            if field_3_i == [number, guy_robin_actual]:
                                y += 1
                        if y == 0:
                            field_sequence_3.append([number, guy_robin_actual])
                            num_cycle += 1                
                    # Condition for the next index
                    if guy_robin_actual > guy_robin_temp:
                        guy_robin_temp = guy_robin_actual
                        field_temp[i_1] += 1
                        i_2 = 0
                        if i_1 == 0:
                            a = field_temp[0]
                    else:
                        if i_2 == 0:
                            field_temp[i_1] -= 1           
                        if i_1 < len(field_temp) -1:  # Prevent index out of range
                            if i_2 == 1:
                                field_temp[i_1] -= 1
                            i_1 += 1
                            field_temp[i_1] += 1
                            i_2 += 1
                            if i_2 == 2:
                                field_temp[i_1] -= 1
                        else:
                            break                                    
                # Next step
                index = 2
                # Defining new field_temp_min
                field_next = []
                for i in range(0, len(field_temp)):
                    if field_temp[i] - 1 <= 1:
                        field_next.append(1)
                    else:
                        field_next.append(field_temp[i] - 1)
                field_next.append(1)                            
            if index == 2:
                # Cycle 2
                i_1 = 0
                i_2 = 0
                field_temp = list(field_temp_min)
                guy_robin_temp = 0.2
                field_temp[0] = a - 1            
                while i_2 < 2 and num_cycle < count_n:            
                    k = 0
                    number = 1
                    decomposition = []
                    for field_i in field_temp:
                        number *= field_primes_init[k] ** field_i
                        k += 1                    
                    # Call decomposition function    
                    decomposition = decomp_number(number)
                    # Call new sigma function, Robopol's theorem
                    sigma_value = sigma_new(decomposition)
                    # Calculate Guy Robin index for sequence
                    if sigma_value > 10 ** 300:
                        logarithm = math.log(math.log(number))       
                        divis = multipl(logarithm, number)                
                        guy_robin_actual = longDivision(sigma_value, divis)                
                    else:
                        guy_robin_actual = sigma_value / (number * math.log(math.log(number)))                    
                    if guy_robin_actual >= 1.72:
                        y = 0
                        for field_3_i in field_sequence_3:
                            if field_3_i == [number, guy_robin_actual]:
                                y += 1
                        if y == 0:
                            field_sequence_3.append([number, guy_robin_actual])
                            num_cycle += 1
                    # Condition for the next index
                    if guy_robin_actual > guy_robin_temp:
                        guy_robin_temp = guy_robin_actual
                        field_temp[i_1] += 1
                        i_2 = 0
                    else:
                        if i_1 > 0:
                            if i_2 == 0:
                                field_temp[i_1] -= 1           
                            if i_1 < len(field_temp) -1:  # Prevent index out of range
                                if i_2 == 1:
                                    field_temp[i_1] -= 1
                                i_1 += 1
                                field_temp[i_1] += 1
                                i_2 += 1
                                if i_2 == 2:
                                    field_temp[i_1] -= 1
                            else:
                                break                                                        
                        else:
                            i_1 = 1
                            field_temp[i_1] += 1
                            guy_robin_temp = 0.2                        
                # Next step
                index = 3                
                # Initialize field_temp_min
                field_temp_min = list(field_next)              
    # Arrange field_sequence_3
    field_sequence_3.sort()
    return field_sequence_3

# Function to update the GUI with outputs
def update_output(field_sequence, sequence_number):
    # Clear previous outputs
    for widget in output_frame.winfo_children():
        widget.destroy()
    
    # Display text output in a Text widget
    text_output = Text(output_frame, height=10, bg='#2b2b2b', fg='white', wrap=WORD)
    text_output.pack(fill=BOTH, expand=True)
    if sequence_number == 1:
        text_output.insert(END, "Field sequence '1' is:\n")
    elif sequence_number == 2:
        text_output.insert(END, "Field sequence '2' is:\n")
    elif sequence_number == 3:
        text_output.insert(END, "Field sequence '3' is:\n")
    # Convert list to string for display with index
    sequence_str = "\n".join([f"{idx+1}: {item[1]:.6f}" for idx, item in enumerate(field_sequence)])
    text_output.insert(END, sequence_str)
    text_output.config(state=DISABLED)
    
    # Plot the sequence graph
    field_chart = [item[-1] for item in field_sequence]
    fig = plt.Figure(figsize=(8, 5), dpi=100)
    ax = fig.add_subplot(111)
    if sequence_number == 1:
        ax.set_title("Guy Robin Index for Sequence '1'", color='white')
    elif sequence_number == 2:
        ax.set_title("Guy Robin Index for Sequence '2'", color='white')
    elif sequence_number == 3:
        ax.set_title("Guy Robin Index for Sequence '3' - Ideal Numbers", color='white')
    ax.grid(True, color='gray')
    ax.set_xlabel("Index in Sequence", color='white')
    ax.set_ylabel("Guy Robin Index Value", color='white')
    ax.plot(field_chart, color='cyan')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    fig.patch.set_facecolor('#2b2b2b')
    ax.set_facecolor('#2b2b2b')

    # Embed the plot in the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=output_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=BOTH, expand=True)

# Function to handle sequence calculation in a separate thread
def run_sequence():
    count_n = get_input()
    if count_n is None:
        return
    selected_sequence = r.get()
    if selected_sequence not in [1, 2, 3]:
        messagebox.showwarning(title="Warning", message="Please select a sequence.")
        return

    # Show progress indicator
    progress_label.config(text="Calculations in progress...")
    start_button.config(state=DISABLED)
    
    # Function to perform calculation and update GUI
    def calculate():
        field_sequence = []
        # clear previous outputs
        for widget in output_frame.winfo_children():
            widget.destroy()
        if selected_sequence == 1:
            field_sequence = sequence_1(count_n)
        elif selected_sequence == 2:
            field_sequence = sequence_2(count_n)
        elif selected_sequence == 3:
            field_sequence = sequence_3(count_n)
        
        # Update the GUI with outputs
        root.after(0, update_output, field_sequence, selected_sequence)
        
        # Hide progress indicator
        root.after(0, lambda: progress_label.config(text=""))
        root.after(0, lambda: start_button.config(state=NORMAL))
    
    # Start calculation in a separate thread
    threading.Thread(target=calculate).start()

# Main loop
root.mainloop()
