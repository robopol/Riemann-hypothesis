import math
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

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

# Defining the necessary constants.
e_gama = 1.7909733665348811333619013505910

# Function to divide large numbers
def longDivision(number, divisor):
    ans = str(number // divisor)
    remainder = number % divisor
    if remainder != 0:
        ans += "."
        for _ in range(10):
            remainder *= 10
            ans += str(remainder // divisor)
            remainder = remainder % divisor
    return float(ans)

# Function: multiplication large numbers
def multipl(float_num, number):
    summ_temp = 0
    point = 0
    int_numb = ""
    divis = 0
    temp = str(float_num)
    for i in range(len(temp)):
        if temp[i] == ".":
            point = i
            break
        else:
            int_numb += temp[i]
    for j in range(point + 1, len(temp)):
        k = (int(temp[j]) * int(number)) // 10 ** (j - point)
        summ_temp += int(k)
    if int(temp[0]) == 0:
        divis = summ_temp
    if int(temp[0]) > 0:
        divis = number * int(int_numb) + summ_temp
    return divis

# Function basic sequence primes
def basic_sequence(n, collect_intermediate=False):
    # Defining the necessary constants   
    k = 1
    sigma_m = 2
    cislo = 6
    index = 2
    field_primes_count_n = [2, 3]
    last_prime = 3
    intermediate_values = []

    target_steps = 10
    step_interval = max(1, n // target_steps)
    current_step = step_interval
    steps_collected = 0

    while index < n and len(intermediate_values) < target_steps:
        divis_1 = 6 * k - 1
        divis_2 = 6 * k + 1
        vysledok1 = 1
        vysledok2 = 1
        hranica = round(math.sqrt(divis_2) + 1)        

        for field in field_primes_count_n:
            if divis_1 % field == 0:
                vysledok1 = 0
            if divis_2 % field == 0:
                vysledok2 = 0    
            if field > hranica:
                break                

        if vysledok1 == 1:
            last_prime = divis_1
            field_primes_count_n.append(divis_1)
            cislo *= divis_1
            sigma_m *= (divis_1 - 1)
            index += 1
            steps_collected += 1
            if collect_intermediate and index >= current_step:
                ratio = cislo / sigma_m
                intermediate_values.append({
                    'last_prime': divis_1,
                    'log_lastprime': math.log(divis_1),
                    'sigma_max_n_over_N': ratio,
                    'e_gamma_log_lastprime': e_gama * math.log(divis_1)
                })
                current_step += step_interval
                if len(intermediate_values) >= target_steps:
                    break

        if vysledok2 == 1 and index < n and len(intermediate_values) < target_steps:
            last_prime = divis_2
            field_primes_count_n.append(divis_2)
            cislo *= divis_2
            sigma_m *= (divis_2 - 1)
            index += 1
            steps_collected += 1
            if collect_intermediate and index >= current_step:
                ratio = cislo / sigma_m
                intermediate_values.append({
                    'last_prime': divis_2,
                    'log_lastprime': math.log(divis_2),
                    'sigma_max_n_over_N': ratio,
                    'e_gamma_log_lastprime': e_gama * math.log(divis_2)
                })
                current_step += step_interval
                if len(intermediate_values) >= target_steps:
                    break
        k += 1                

    # Calculate final values
    log_lastprime = math.log(last_prime)
    if cislo >= 10**300:
        ratio = longDivision(cislo, sigma_m)
    else:
        ratio = cislo / sigma_m

    final_values = {
        'last_prime': last_prime,
        'log_lastprime': log_lastprime,
        'sigma_max_n_over_N': ratio,
        'e_gamma_log_lastprime': e_gama * log_lastprime
    }

    return final_values, intermediate_values

# Tkinter Application
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("σ_max(N) for Highly Composite Numbers")
        self.geometry("1000x800")
        self.minsize(800, 600)
        
        # Flag to indicate if the window is closing
        self.closing = False

        # Configure grid layout
        self.columnconfigure(0, weight=1)
        self.rowconfigure(2, weight=1)

        # Define styles
        style = ttk.Style(self)
        style.configure("Treeview", font=("Arial", 10))
        style.configure("Treeview.Heading", font=("Arial", 12, "bold"))
        style.configure("TLabel", font=("Arial", 12))
        style.configure("TButton", font=("Arial", 12))
        style.configure("TEntry", font=("Arial", 12))

        # Input Frame
        input_frame = ttk.Frame(self)
        input_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        input_frame.columnconfigure(1, weight=1)

        ttk.Label(input_frame, text="Enter number of primes (n):").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.n_entry = ttk.Entry(input_frame)
        self.n_entry.grid(row=0, column=1, padx=5, pady=5, sticky="ew")
        self.n_entry.focus()

        compute_button = ttk.Button(input_frame, text="Start", command=self.start_computation)
        compute_button.grid(row=0, column=2, padx=5, pady=5, sticky="e")

        # Results Frame
        results_frame = ttk.LabelFrame(self, text="Results")
        results_frame.grid(row=1, column=0, padx=10, pady=10, sticky="ew")

        columns = ("step", "last_prime", "log_lastprime", "sigma_max_n_over_N", "e_gamma_log_lastprime")
        self.tree = ttk.Treeview(results_frame, columns=columns, show="headings", height=10)
        self.tree.heading("step", text="Step")
        self.tree.heading("last_prime", text="Last Prime")
        self.tree.heading("log_lastprime", text="log(Last Prime)")
        self.tree.heading("sigma_max_n_over_N", text="σ_max(n)/N")
        self.tree.heading("e_gamma_log_lastprime", text="e^γ * log(Last Prime)")

        self.tree.column("step", width=60, anchor="center")
        self.tree.column("last_prime", width=150, anchor="center")
        self.tree.column("log_lastprime", width=200, anchor="center")
        self.tree.column("sigma_max_n_over_N", width=150, anchor="center")
        self.tree.column("e_gamma_log_lastprime", width=200, anchor="center")

        # Adding scrollbar to the treeview
        scrollbar = ttk.Scrollbar(results_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        self.tree.pack(side="left", fill="both", expand=True)

        # Graph Frame
        graph_frame = ttk.LabelFrame(self, text="Graph of σ_max(n)/N and e^γ * log(Last Prime)")
        graph_frame.grid(row=2, column=0, padx=10, pady=10, sticky="nsew")
        graph_frame.columnconfigure(0, weight=1)
        graph_frame.rowconfigure(0, weight=1)
        
        self.figure, self.ax = plt.subplots(figsize=(8,4))
        self.canvas = FigureCanvasTkAgg(self.figure, master=graph_frame)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        # Bind the window close event
        self.protocol("WM_DELETE_WINDOW", self.on_closing)

    def start_computation(self):
        if self.closing:
            return  # Prevent starting computation if closing

        n_str = self.n_entry.get()
        try:
            n = int(n_str)
            if n <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid integer greater than 0.")
            return

        # Clear previous results
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.ax.clear()
        self.canvas.draw()

        # Disable the Start button to prevent multiple clicks
        compute_button = self.children['!frame'].children['!button']
        compute_button.config(state='disabled')

        try:
            # Perform the computation
            final_values, intermediate = basic_sequence(n, collect_intermediate=True)
            
            # Select the last 10 intermediate values
            sampled = intermediate[:10]  # Since we collect at step_interval, it's already 10

            # Prepare data for graph
            step_size = n // 10
            x_values = [step_size * i for i in range(1, len(sampled)+1)]
            y_values_sigma = [item['sigma_max_n_over_N'] for item in sampled]
            y_values_e_gamma = [item['e_gamma_log_lastprime'] for item in sampled]

            # Insert into table
            for idx, data in enumerate(sampled, start=1):
                self.tree.insert("", "end", values=(
                    idx,
                    data['last_prime'],
                    f"{data['log_lastprime']:.10f}",
                    f"{data['sigma_max_n_over_N']:.10f}",
                    f"{data['e_gamma_log_lastprime']:.10f}"
                ))

            # Plot σ_max(n)/N
            self.ax.plot(x_values, y_values_sigma, marker='o', linestyle='-', color='blue', label='σ_max(n)/N')
            # Plot e^γ * log(Last Prime)
            self.ax.plot(x_values, y_values_e_gamma, marker='s', linestyle='--', color='red', label='e^γ * log(Last Prime)')

            self.ax.set_title("σ_max(n)/N and e^γ * log(Last Prime) vs n")
            self.ax.set_xlabel("n")
            self.ax.set_ylabel("Values")
            self.ax.legend()
            self.ax.grid(True)
            self.figure.tight_layout()
            self.canvas.draw()

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during computation:\n{e}")
        finally:
            # Re-enable the Start button
            compute_button.config(state='normal')

    def on_closing(self):
        self.closing = True
        self.destroy()

# Run the application
if __name__ == "__main__":
    app = App()
    app.mainloop()
