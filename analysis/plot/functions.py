""" head functions.py -n 2
This is a script of functions that I couldn't get to work from a different package. It's not meant to be run
"""

# %% Import packages
import math


# %% Import working test
def test():
    print("working")


# %% Chi Squared Calculations
def chi_squared(observed_tab, expected_tab = None, test_type = None):
    """
    Calculates the chi-squared statistic and corresponding p-value for the input table

    Parameters
    ----------
    observed_tab : list of lists
        The observed values in table format.
    expected_tab : list of lists, optional
        The expected values in table format. Must be the same shape as observed_tab if it exists. The default is None.
    test_type : string, optional
        If it is the independence test, put 'independence'. The default is even chance of all options.

    Returns
    -------
    chi_squared_statistic : float
        The chi-squared value from the chi-squared test.
    p_value : float
        The p-value from the chi-squared test.

    """
    
    # Rearrange the observed variable from a table into a vector
    observed = [] # Initialize the vector
    for row in observed_tab: # for each row in the table
        for entry in row: # for each entry in the row
            observed.append(entry) # put that entry into the vector
    
    # Find the total and df
    total_value = sum(observed) # the total of the table is the sum of the vector
    df = (len(observed_tab) - 1) * (len(observed_tab[0]) - 1) # the degrees of freedom is (number of rows-1) * (number of columns-1)
    
    if expected_tab == None and test_type == None: # If there is no listed expected value and no listed testing type
        # Assume that you're comparing it to something where there is even chance for each entry in the table
        expected = [] # Start a vector for the expected values
        for ind in range(len(observed)): # for each entry in the observed vector
            expected.append(total_value/len(observed)) # add an expected value that is the average of the observed values
    
    elif expected_tab == None and test_type == "independence": # If the test type is specified as 'independent'
        # The expected value is dependent on the percentage of the total that is in each row and column overall
        
        row_totals = [] # start a new vector for the total value in each row
        column_totals = [0] * len(observed_tab[0]) # start a new vector for the total value of each column
        
        for row in observed_tab: # Look through each row in the observed table
            row_totals.append(sum(row)) # the total of the row is the sum of the row
            for col_ind in range(len(row)): # loop through the values in the current row
                column_totals[col_ind] += row[col_ind] # add the current value to the corersponding column total
                
        expected_tab = [] # Start a table for the expected values
        for row_val in row_totals: # loop through each row in the total values
            expected_tab.append([]) # start a new row in the  expected for the row from the observed
            for col_val in column_totals: # for each value in the columns
                exp_temp = col_val * row_val / total_value # caluclate the expected value for that row/column pair
                expected_tab[-1].append(exp_temp) # Put that calculated value into the table
        
        # Rearrange the expected variable
        expected = [] # Start a vector for the expected table
        for row in expected_tab: # For each row in the table...
            for entry in row: # ...and each entry in that row...
                expected.append(entry) # ...put the entry into the vector
    
    component = [] # the component vector holds the chi-squared value for each entry in the original table, matching the vector orders of the variables expected and observed
    res_all = [[],[]] # the residuals and the squares of the residuals stored from their temporary variables
    
    # Loop through each observed value
    for index in range(len(observed)):
        res = observed[index] - expected[index] # the residual is the difference of the observed and expected value
        res_all[0].append(res) # save the residual from the temporary variable
        res2 = res * res # the square of the residual
        res_all[1].append(res2) # save the square of the residual from the temporary variable
        com = res2 / expected[index] # the chi-squared component is the square of residuals divided by the expected
        component.append(com) # save the component of the chi-squared statistic for the current pair of expected and observed values
    
    chi_squared_statistic = sum(component) # the total chi-squared statistic is the sum of the components for each observed value
    
    # Ideally the CDF would be defined here and then the p-value would be calculated
    # Sadly, this is not an ideal world
    
    return chi_squared_statistic#, p_value


# %% Linear model fitting
def Linear_Fitting(X, Y):
    """
    Caluclates the slope and intercept of a best fit line with the given X and Y vectors using least-squares regression

    Parameters
    ----------
    X : list
        The x-values to be fit to a linear regression model.
    Y : list
        The y-values to be fit to a linear regression model.

    Returns
    -------
    m : float
        The slope of the best-fit line.
    b : float
        The intercept of the best-fit line.

    """
    
    X2 = [] # initialize a list for each X^2 value
    XY = [] # initialize a list for each X*Y value
    
    for ind in range(len(X)): # Loop through each x and y value within the inptu vectors
        X2.append(X[ind] * X[ind]) # add the x^2 value to the correct vector
        XY.append(X[ind] * Y[ind]) # add the x*y value to the correct vector
    
    # Find the sums of each of these values calcualted and the total number of observations
    X_s = sum(X)
    Y_s = sum(Y)
    X2_s = sum(X2)
    XY_s = sum(XY)
    n = len(X)
    
    # These are the equations for linear regression using the least squares method
    # these were taken from Khan Academy
    m_num = (n * XY_s) - (X_s * Y_s)
    m_denom = (n * X2_s) - (X_s * X_s)
    m = m_num / m_denom
    
    b = (Y_s - (m * X_s) ) / n
    
    # return the slope and itnercept
    return m, b
    

# %% Linear model assessment
def R_squared_calculation(Y,Y_predicted):
    """
    

    Parameters
    ----------
    Y : list
        The observed Y values from a linear regression model.
    Y_predicted : list
        The predicted Y values from a linear regression model.

    Returns
    -------
    R_squared : float
        The coefficient of determination using the original and predicted y values.

    """
    
    
    #I don't know where the equation for this came from
    mean_y = sum(Y) / len(Y) # find the mean of the y values
    
    # initialize the numerator and denominator of the R_squared value
    SSR = 0
    SST = 0
    SSR_all = []
    SST_all = []
    
    # Calculate the sum of squares of residuals and the total sum of squares
    for index in range(len(Y)):
        ssr_temp = Y_predicted[index] - Y[index]
        sst_temp = Y[index] - mean_y
        
        SSR_all.append(ssr_temp)
        SST_all.append(sst_temp)
        
        SSR += ssr_temp * ssr_temp
        SST += sst_temp * sst_temp
        
    # Use the sums of squares to calculate the R squareed value
    R_squared = 1 - (SSR / SST)
    
    
    """
    Nope
    
    x2 = []
    y2 = []
    xy = []
    n = len(Y)
    
    for index in range(len(Y)):
        x2.append(Y[index]*Y[index])
        y2.append(Y_predicted[index]*Y_predicted[index])
        xy.append(Y[index]*Y_predicted[index])
    
    x2_sum = sum(x2)
    y2_sum = sum(y2)
    x_sum = sum(Y)
    y_sum = sum(Y_predicted)
    xy_sum = sum(xy)
    
    numerator = (n * xy_sum) - (x_sum * y_sum)
    den1 = math.sqrt( (n*x2_sum) - (x_sum * x_sum) )
    den2 = math.sqrt( (n*y2_sum) - (y_sum * y_sum) )
    denominator = den1 + den2
    
    Correlation_Coefficient = numerator/denominator
    R_squared = Correlation_Coefficient * Correlation_Coefficient
    """
    """
    numerator = []
    Y_avg = sum(Y)/len(Y)
    Y_pred_avg = sum(Y_predicted)/len(Y_predicted)
    
    for index in range(len(Y)):
        numerator.append( (Y[index]-Y_avg) * (Y_predicted[index]-Y_pred_avg))
    denominator = ( (len(Y)-1) * () * ())
    """
    return R_squared




def squared_error(ys_orig,ys_line):
    answer = 0
    for i in range(len(ys_orig)):
        answer += (ys_line[i] - ys_orig[i]) * (ys_line[i] - ys_orig[i])
    return answer

def coefficient_of_determination(ys_orig,ys_line):
    y_mean_line = [sum(ys_orig)/len(ys_orig) for y in ys_orig]
    squared_error_regr = squared_error(ys_orig, ys_line)
    squared_error_y_mean = squared_error(ys_orig, y_mean_line)
    return 1 - (squared_error_regr/squared_error_y_mean)
