from PIL import Image
 
def convertImage():
    img = Image.open("/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Testing/FTD_GM.png")
    img = img.convert("RGBA")
 
    datas = img.getdata()
 
    newData = []
 
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            newData.append((255, 255, 255, 0))
        else:
            newData.append(item)
 
    img.putdata(newData)
    img.save("/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Output_Python_NewFTD/New.png", "PNG")
    print("Successful")
 
convertImage()