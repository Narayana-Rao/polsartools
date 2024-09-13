import matplotlib.pyplot as plt

def plot_data(data_array):
    """
    Visualizes the data array using matplotlib.
    """
    plt.imshow(data_array, cmap='gray')
    plt.colorbar()
    plt.show()
