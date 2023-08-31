import random, argparse

# IMPORTANT: set the n_pts arg in this file before building -> vorocrust-meshing/tests/datasets/OneBox_monitoring_points/run_series.sh (line number 4), if you want to use this with ctest
# To use help menu -> $ python3 check_points_in_mesh.py -h
def arg_setup():
	parser = argparse.ArgumentParser()
	parser.add_argument("--num_pts", action="store", dest="num_pts", help="This is the number of points you want to randomly pull from the closest seeds csv file to check for existance in mesh.vcg", type=int)

	return parser.parse_args()

# Specifying exception for if all points don't make it into the final mesh
class NotAllPointsInFinalMeshError(Exception):
	def __init__(self, message="All points in monitoring_points_closest_seeds.csv did not make it into mesh.vcg, test FAIL"):
		self.message = message
		super().__init__(self.message)

# This will validate that we are not specifying a number of points that is greater than the total number in the file
# If user enters a value greater than the number of points in the file, it will be set to the max number of points
def arg_validation_check(num, tnp):
	if (num > tnp):
		h = tnp
	else:
		h = num

	return h

# This function will return a list of random points from the monitoring_points_closest_seeds file
# The number of points is specified by the user
def get_random_points_from_file(n, lines, tnp):
	line_number = []
	points = []
	
	while(len(line_number) < n):
		temp_rand = random.randint(1, tnp) # We need to start at 1 since the 0 line is header info
		if temp_rand not in line_number:
			line_number.append(temp_rand)

	for each in line_number:
		points.append(lines[each])

	return points

# This function will put the points in a format that can be used to match the points if they
# exist in mesh.vcg, this cleanup includes getting rid of the commas and the radius data
def cleanup_points(points):
	clean = []
	for x in range(0, len(points)):
		temp_split = points[x].split(",")
		coords = temp_split[:3]
		clean.append(coords[0] + coords[1] + coords[2])

	return clean

# This function searches mesh.vcg for the points, if found, a counter is incremented
# returned value is the percent of points found
def search_mesh_vcg(sf, n, cl):
	found = 0
	for x in range(0, n):
		#print(cl[x])
		if cl[x] in sf:
			found = found + 1

	return found / n

# This function will verify the results for us, if not all points made it into the final
# mesh then raise error and the test will fail
def check_percentage(p):
	if (p != 1.0):
		raise NotAllPointsInFinalMeshError
	else:
		return

# Main function
if __name__ == "__main__":
	args = arg_setup()
	mvcg = open("mesh.vcg", "r")
	mpcs = open("monitoring_points_closest_seeds.csv", "r")

	search_file = mvcg.read()
	lines_mpcs = mpcs.readlines()
	total_num_points = len(lines_mpcs) - 1

	how_many_points = arg_validation_check(args.num_pts, total_num_points)
	points = get_random_points_from_file(how_many_points, lines_mpcs, total_num_points)

	cl = cleanup_points(points)

	percent_found = search_mesh_vcg(search_file, how_many_points, cl)

	check_percentage(percent_found)


	
