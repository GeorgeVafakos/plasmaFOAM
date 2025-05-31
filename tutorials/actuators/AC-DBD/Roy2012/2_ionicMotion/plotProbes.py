import matplotlib.pyplot as plt

# Function to read the data file
def read_and_plot(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Skip comment lines and extract numerical data
    data_lines = [line.strip() for line in lines if not line.startswith('#')]

    # Convert the extracted lines into a 2D list of floats
    data = [list(map(float, line.split())) for line in data_lines]

    # Separate the first column (x-axis) and the rest (y-axis)
    x_values = [row[0] for row in data]
    y_values = [row[1:] for row in data]

    # Transpose y_values to separate each column
    y_columns = list(zip(*y_values))

    # Plot each y column against the x column
    for i, y in enumerate(y_columns):
        plt.plot(x_values, y, '-', label=f'Column {i+1}')

    plt.xlabel('Time (s)')
    plt.ylabel('rhoq (C/m3)')
    plt.title('Charge density at x=0.1mm and y=0.1mm')
    plt.legend()
    plt.grid()
    plt.show()

# Provide the file path
file_path = './postProcessing/probes/0/rhoq'
read_and_plot(file_path)

