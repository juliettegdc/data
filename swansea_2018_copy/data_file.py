# Import necessary modules
import geopandas as gpd
import shapely.geometry
import utm

# Read file using gpd.read_file()
data_GBS = gpd.read_file("/data/BGS/DigSBS/Data/BGS_250k_SeaBedSediments_WGS84_v3.shp")


outfile = open('MODELandBGS.txt', 'w')
outfile.write('x y lat lon ModelMeanBSS ModelMaxBSS BGS-Folk \n')

file_mean = open('data_binned_shear.csv', 'r')
file_max = open("data_binned_max_shear.csv", 'r')

lines_max=file_max.readlines()
#print lines[25]
#file_mean.readline() #jump first line
#file_max.readline()

polygon = shapely.geometry.Polygon([(-5.160, 51.876), (-2.353, 51.876), (-5.160, 50.526), (-2.353, 50.526)])

#for _ in range(1980):
      #  next(file_mean)
cnt = 0
for line in file_mean:

    #add error message if files don't have same length

    b = line.split('\t')
    c = lines_max[cnt].split('\t')

    x = float(b[0])
    y = float(b[1])
    number_mean = float(b[2])
    number_max = float(c[2])

    lat, lon = utm.to_latlon(x, y, 30, 'U')  # converts x,y of mesh into lon, lat

    point = shapely.geometry.Point(lon, lat)
    print(point, (cnt+1))

    cnt += 1
    check = 0

    for index, row in data_GBS.iterrows():  # selection.iterrows():
        # print(index, row['geometry'])

        shape = shapely.geometry.asShape(row['geometry'])

        #if shape.within(polygon):
        if point.within(shape):
            print("Found shape for point.")
            print(index, row['RCS'])  # row['geometry'],
            outfile.write(
                str(x) + '\t' + str(y) + '\t' + str(lat) + '\t' + str(lon) + '\t' + str(number_mean) + '\t' + str(
                    number_max) + '\t' + row['RCS'] + '\n')
            check = 1
            break

    if check == 0:
        outfile.write(
            str(x) + '\t' + str(y) + '\t' + str(lat) + '\t' + str(lon) + '\t' + str(number_mean) + '\t' + str(
                number_max) + '\t' + 'none ' + '\n')


outfile.close()
file_mean.close()
file_max.close()

