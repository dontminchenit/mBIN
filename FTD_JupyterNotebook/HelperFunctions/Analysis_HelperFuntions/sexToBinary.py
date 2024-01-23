def sexLUT(sexstr):
    """convert string (gender type) to int

    Args:
        sexstr (str): Male or Female

    Raises:
        Exception: Other than Male or Female

    Returns:
        0 or 1 (int): 0 - Male / 1 - Female
    """
    if sexstr == 'Male':
        return 0
    elif sexstr == 'Female':
        return 1
    else:
        raise Exception("Problem with converting sex to integer")

def sexToBinary(listOfValues):
    return [sexLUT(i) for i in listOfValues]





def ADNCLUT(ADNCstr):
    """convert string (ADNC) to int

    Args:
        ADNCstr (str): a_Not, b_Low, c_Int, d_High

    Raises:
        Exception: Other than a_Not, b_Low, c_Int, d_High

    Returns:
        0 or 1 (int): 0 - a_Not / 1 - b_Low / 2 - c_Int / 3 - d_High
    """
    if ADNCstr == 'a_Not':
        return 0
    elif ADNCstr == 'b_Low':
        return 1
    elif ADNCstr == 'c_Int':
        return 2
    elif ADNCstr == 'd_High':
        return 3
    else:
        raise Exception("Problem with converting sex to integer")

def ADNCToBinary(listOfValues):
    return [ADNCLUT(i) for i in listOfValues]