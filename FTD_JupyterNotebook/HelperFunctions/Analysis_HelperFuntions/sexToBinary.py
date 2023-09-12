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