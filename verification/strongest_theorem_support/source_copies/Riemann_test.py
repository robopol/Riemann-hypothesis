# Importing libraries
import math
import threading
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')  # Setting the backend for Tkinter

# Function for dividing large numbers
def longDivision(number, divisor):
    ans = str(number // divisor)
    remainder = 1
    if number % divisor != 0:
        ans += "."
        remainder = number % divisor
    for _ in range(1, 11):
        temp = remainder * 10
        ans += str(temp // divisor)
        remainder = temp % divisor
    try:
        ans = float(ans)
    except:
        ans = 0.0
    return ans

# Function for multiplying large numbers
def multipl(float_num, number):
    summ_temp = 0
    point = 0
    int_numb = ""
    multiple_sum = 0
    temp = str(float_num)
    for i in range(0, len(temp) - 1):
        if temp[i] == ".":
            point = i
            break
        else:
            int_numb += temp[i]
    for j in range(point + 1, len(temp)):
        k = (int(temp[j]) * int(number)) // (10 ** (j - point))
        summ_temp += int(k)
    if int(temp[0]) == 0:
        multiple_sum = summ_temp
    elif int(temp[0]) > 0:
        multiple_sum = number * int(int_numb) + summ_temp
    return multiple_sum

# Function for calculating the sequence of ideal numbers
def basic_sequence(n, step_size):
    # Defining necessary constants
    k = 1
    sigma_seq_1 = 12
    numb_seq_1 = 6
    index = 2
    field_primes_count_n = [2, 3]
    last_prime = 3  # Initialization of the last prime number
    ratio = sigma_seq_1 / numb_seq_1
    num_max = numb_seq_1

    while index < n:
        divis_1 = 6 * k - 1
        divis_2 = 6 * k + 1
        result1 = 1
        result2 = 1
        limit = round(math.sqrt(divis_2) + 1)
        for field in field_primes_count_n:
            if divis_1 % field == 0:
                result1 = 0
            if divis_2 % field == 0:
                result2 = 0
            if field > limit:
                break
        if result1 == 1:
            last_prime = divis_1
            field_primes_count_n.append(divis_1)
            numb_seq_1 *= divis_1
            sigma_seq_1 *= (divis_1 + 1)
            index += 1
        if result2 == 1:
            last_prime = divis_2
            field_primes_count_n.append(divis_2)
            numb_seq_1 *= divis_2
            sigma_seq_1 *= (divis_2 + 1)
            index += 1
        k += 1

    # Calculation of the Guy Robin index
    if sigma_seq_1 >= 10**300:
        logarithm = math.log(math.log(numb_seq_1))
        divisor = multipl(logarithm, numb_seq_1)
        guy_robin = longDivision(sigma_seq_1, divisor)
    else:
        guy_robin = sigma_seq_1 / (numb_seq_1 * math.log(math.log(numb_seq_1)))

    # Calculation for the ideal number
    constant_1 = 0
    constant_2 = 0
    guy_max = guy_robin
    seq_number = 0
    number = 2
    index_1 = 0
    sigma_new = sigma_seq_1
    sigma_aux = 1
    numb_seq_3 = numb_seq_1

    while constant_2 < 2 and index_1 < len(field_primes_count_n):
        sigma_aux = 1
        for field in field_primes_count_n:
            if field != number:
                sigma_aux *= (field + 1)
        power = 1
        auxiliary_number = field_primes_count_n[index_1]
        while constant_1 < 1:
            power += 1
            sigma_new += (number ** power) * sigma_aux
            numb_seq_3 *= number
            auxiliary_number += number ** power
            # Calculation of the Guy Robin index
            if sigma_new >= 10**300:
                logarithm = math.log(math.log(numb_seq_3))
                divisor = multipl(logarithm, numb_seq_3)
                guy_robin = longDivision(sigma_new, divisor)
            else:
                guy_robin = sigma_new / (numb_seq_3 * math.log(math.log(numb_seq_3)))
            # Conditions for guy_max
            if guy_robin > guy_max:
                constant_1 = 0
                constant_2 = 0
                guy_max = guy_robin
                num_max = numb_seq_3
                ratio = longDivision(sigma_new, num_max)
            else:
                constant_1 = 1
                constant_2 += 1
        # Change of field
        field_primes_count_n[index_1] = auxiliary_number
        seq_number += 1
        index_1 += 1
        if index_1 >= len(field_primes_count_n):
            break  # Preventing IndexError
        number = field_primes_count_n[index_1]
        constant_1 = 0

    # Output conditions
    log_lastprime = math.log(last_prime)
    logN = math.log(num_max)
    # Returning the results instead of printing them
    return guy_max, log_lastprime, ratio, logN, num_max

# Definition of the GUI application
class RiemannHypothesisApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Riemann Hypothesis - advanced testing")
        self.root.geometry("1000x1000")
        self.root.configure(bg='#2b2b2b')  # Dark background

        # Configuring styles
        self.style = ttk.Style()
        self.style.theme_use('clam')
        self.style.configure('.', background='#2b2b2b', foreground='white', fieldbackground='#444444')
        self.style.configure('TLabel', background='#2b2b2b', foreground='white')
        self.style.configure('TButton', background='#444444', foreground='white')
        self.style.configure('TEntry', fieldbackground='#444444', foreground='white')
        self.style.configure('TFrame', background='#2b2b2b')

        # Creating widgets
        self.create_widgets()

    def create_widgets(self):
        # Frame for inputs
        input_frame = ttk.LabelFrame(self.root, text="Input Parameters", padding=(20, 20))
        input_frame.pack(padx=20, pady=10, fill=tk.X)

        # Input for the maximum number
        ttk.Label(input_frame, text="Enter the maximum number (e.g., 10000):").pack(anchor=tk.W, pady=(0, 5))
        self.num_entry = ttk.Entry(input_frame, width=20)
        self.num_entry.pack(anchor=tk.W, pady=(0, 10))
        self.num_entry.insert(0, "10000")  # Default value

        # Button to start the computation
        self.start_button = ttk.Button(input_frame, text="Start Computation", command=self.start_computation)
        self.start_button.pack(pady=10)

        # Frame for results
        results_frame = ttk.LabelFrame(self.root, text="Results", padding=(20, 20))
        results_frame.pack(padx=20, pady=10, fill=tk.BOTH, expand=True)

        # Text box for displaying results
        self.results_text = tk.Text(results_frame, bg='#2b2b2b', fg='white', wrap=tk.WORD, height=15)
        self.results_text.pack(fill=tk.BOTH, expand=True)
        self.results_text.config(state=tk.DISABLED)

        # Frame for graph
        plot_frame = ttk.LabelFrame(self.root, text="Guy Robin Index Graph", padding=(20, 20))
        plot_frame.pack(padx=20, pady=10, fill=tk.BOTH, expand=True)

        # Matplotlib Figure
        self.fig, self.ax = plt.subplots(figsize=(8, 4), dpi=100)
        self.ax.set_title("Guy Robin Index Progress", color='white')
        self.ax.set_xlabel("Computation Step", color='white')
        self.ax.set_ylabel("Guy Robin Index", color='white')
        self.ax.grid(True, color='gray')
        self.ax.set_facecolor('#2b2b2b')
        self.fig.patch.set_facecolor('#2b2b2b')
        self.line, = self.ax.plot([], [], color='cyan')

        # Embedding the graph in Tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def start_computation(self):
        # Getting input from the user
        try:
            max_number = int(self.num_entry.get())
            if max_number <= 0:
                raise ValueError
        except ValueError:
            messagebox.showwarning("Invalid Input", "Please enter a positive integer for the maximum number.")
            return

        # Disable the button during computation
        self.start_button.config(state=tk.DISABLED)

        # Reset results and graph
        self.results_text.config(state=tk.NORMAL)
        self.results_text.delete(1.0, tk.END)
        self.results_text.config(state=tk.DISABLED)
        self.ax.cla()  # Clearing previous graph
        self.ax.set_title("Guy Robin Index Progress", color='white')
        self.ax.set_xlabel("Computation Step", color='white')
        self.ax.set_ylabel("Guy Robin Index", color='white')
        self.ax.grid(True, color='gray')
        self.ax.set_facecolor('#2b2b2b')
        self.fig.patch.set_facecolor('#2b2b2b')
        self.line, = self.ax.plot([], [], color='cyan')
        self.canvas.draw()

        # Starting computation in a separate thread
        threading.Thread(target=self.compute_sequence, args=(max_number,), daemon=True).start()

    def compute_sequence(self, max_number):
        results = []
        x_vals = []
        y_vals = []

        # Defining the step size
        step_size = max(1, math.ceil(max_number / 10))
        total_steps = math.ceil(max_number / step_size)

        for step in range(1, total_steps + 1):
            current_n = step * step_size
            if current_n > max_number:
                current_n = max_number

            try:
                # Performing the computation
                guy_max, log_lastprime, ratio, logN, num_max = basic_sequence(current_n, step_size)

                # Storing the results
                results.append(guy_max)
                x_vals.append(current_n)
                y_vals.append(guy_max)

                # Updating the results text box
                self.results_text.config(state=tk.NORMAL)
                self.results_text.insert(tk.END, f"Step {step} ({current_n}): {guy_max:.6f}\n")
                self.results_text.see(tk.END)
                self.results_text.config(state=tk.DISABLED)

                # Updating the graph
                self.line.set_data(x_vals, y_vals)
                self.ax.relim()
                self.ax.autoscale_view()
                self.canvas.draw()

            except Exception as e:
                messagebox.showerror("Error", f"An error occurred during the computation: {e}")
                break  # Ending the loop on error

        # Re-enable the button after computation
        self.start_button.config(state=tk.NORMAL)

# Running the application
if __name__ == "__main__":
    root = tk.Tk()
    app = RiemannHypothesisApp(root)
    root.mainloop()
