import platform

from support_functions import EXIT_TRACEBACK

class Paths_class:
	def __init__(self):
		self.gmsh          = ''
		self.DPG_ROOT      = ''
		self.meshes        = ''
		self.cases         = ''
		self.control_files = ''

	def set_paths(self,user):
		""" Set paths based on the current user.  """
		if (user == 'PZwan'):
			if any('Darwin' in string for string in platform.uname()):
				self.gmsh     = '/Applications/Gmsh.app/Contents/MacOS/gmsh'
				self.DPG_ROOT = '/Users/philip/Desktop/DPGSolver/'
			else:
				# Ubuntu (home)
				self.gmsh = '/home/philip/Desktop/research/programs/gmsh/gmsh-2.12.0-Linux/bin/gmsh'
				self.DPG_ROOT = '/home/philip/Desktop/research/codes/DPGSolver/'
			self.meshes = self.DPG_ROOT+'meshes/'
		else:
			print('Add an option for yourself as a user.')
			EXIT_TRACEBACK()

		self.cases         = self.DPG_ROOT+'cases/'
		self.control_files = self.cases+'control_files/'


def find_MeshType(N,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving):
	Found = 0
	for i in range(0,N):
		if (MeshName.find(MeshCurving[i]+MeshTypes[i]) != -1):
			# Eliminate occurence of finding Curved in ToBeCurved ...
			if (MeshName.find('ToBeCurved') != -1 and MeshCurving[i].find('ToBeCurved') == -1):
				continue

			Found = 1
			MeshType    = [MeshTypes[i]]
			MeshCurving = [MeshCurving[i]]
			if (len(MeshTypesPrefix) == 1):
				MeshTypesPrefix = [MeshTypesPrefix[0]]
			else:
				MeshTypesPrefix = [MeshTypesPrefix[i]]

			break;

	if (Found == 0):
		print("Did not find the MeshType.\n")
		print(MeshName,"\n",MeshCurving,"\n",MeshTypes,"\n")
		EXIT_TRACEBACK()

	return [MeshType, MeshTypesPrefix, MeshCurving]

class TestCase_class:
	""" Class storing relevant data for a control file group."""

	def __init__(self,name):
		self.name    = name
		self.VarName = ''

		self.MeshTypes  = []
		self.Path = ''

		self.GeoDeps     = ''
		self.MeshOutputs = ''

	def add_MeshTypes(self,Paths,MeshName):
		if (self.name.find('update_h') != -1 or
		    self.name.find('L2_proj')  != -1):
			if (self.name.find('update_h') != -1):
				self.VarName = 'UPDATE_H'
			elif (self.name.find('L2_proj_p') != -1):
				self.VarName = 'L2_PROJ_P'
			elif (self.name.find('L2_proj_h') != -1):
				self.VarName = 'L2_PROJ_H'

			NTotal = 6
			MeshTypesPrefix = ['' for i in range(NTotal)]
			MeshCurving     = ['' for i in range(NTotal)]
			MeshTypes       = ['TRI','QUAD','TET','HEX','WEDGE','PYR']
			if   (MeshName.find('all')    != -1):
				iRange = range(0,NTotal)
			else:
				iRange = range(0,1)
				[MeshTypes,MeshTypesPrefix,MeshCurving] = find_MeshType(NTotal,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving)
		elif (self.name.find('linearization') != -1):
			self.VarName = 'LINEARIZATION'

			NTotal = 3
			MeshTypesPrefix = ['' for i in range(NTotal)]
			MeshCurving     = ['ToBeCurved' for i in range(NTotal)]
			MeshTypes       = ['MIXED2D','MIXED3D_TP','MIXED3D_HW']
			if   (MeshName.find('all')    != -1):
				iRange = range(0,NTotal)
			else:
				iRange = range(0,1)
				[MeshTypes,MeshTypesPrefix,MeshCurving] = find_MeshType(NTotal,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving)
		elif (self.name.find('Poisson') != -1):
			if (self.name.lower().find('test') != -1):
				self.VarName = 'POISSON_TEST'

				MeshTypesPrefix = ['n-Ball_HollowSection_','n-Ellipsoid_HollowSection_','n-Ellipsoid_HollowSection_']
				MeshCurving     = ['Curved','Curved','Curved']
				MeshTypes       = ['MIXED2D','TRI','QUAD']
				NTotal = len(MeshTypes)
				if   (MeshName.find('all')    != -1):
					iRange = range(0,NTotal)
				else:
					iRange = range(0,1)
					[MeshTypes,MeshTypesPrefix,MeshCurving] = find_MeshType(NTotal,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving)
			else:
				EXIT_TRACEBACK()
		elif (self.name.find('Euler') != -1):
			if (self.name.lower().find('test') != -1):
				self.VarName = 'EULER_TEST'

				MeshCurving     = ['Curved','ToBeCurved','ToBeCurved','ToBeCurved','ToBeCurved','ToBeCurved','']
				MeshTypes       = ['MIXED2D','MIXED2D','MIXED3D_TP','TET','HEX','WEDGE','QUAD']
				NTotal = len(MeshTypes)
				MeshTypesPrefix = ['SupersonicVortex_' for i in range(0,6)]
				MeshTypesPrefix.append('PeriodicVortex_')
				if   (MeshName.find('all')    != -1):
					iRange = range(0,NTotal)
				else:
					iRange = range(0,1)
					[MeshTypes,MeshTypesPrefix,MeshCurving] = find_MeshType(NTotal,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving)
			else:
				print(self.name)
				EXIT_TRACEBACK()
		elif (self.name.find('NavierStokes') != -1):
			if (self.name.lower().find('test') != -1):
				self.VarName = 'NAVIERSTOKES_TEST'

				MeshCurving     = ['ToBeCurved']
				MeshTypes       = ['TRI']
				NTotal = len(MeshTypes)
				MeshTypesPrefix = ['TaylorCouette_' for i in range(0,NTotal)]
				if   (MeshName.find('all')    != -1):
					iRange = range(0,NTotal)
				else:
					iRange = range(0,1)
					[MeshTypes,MeshTypesPrefix,MeshCurving] = find_MeshType(NTotal,MeshTypes,MeshTypesPrefix,MeshName,MeshCurving)
			else:
				print(self.name)
				EXIT_TRACEBACK()
		else:
			EXIT_TRACEBACK()

		for i in iRange:
			TypeCurrent = MeshType_class(MeshTypes[i],MeshTypesPrefix[i]+MeshCurving[i])

			TypeCurrent.set_parameters(self,Paths)

			self.MeshTypes.append(TypeCurrent)

	def set_paths(self,Paths):
		if (self.name.find('update_h') != -1):
			self.name = 'Test_update_h_'
			self.Path = Paths.control_files+'test/update_h/'
		elif (self.name.find('L2_proj_p') != -1):
			self.name = 'Test_L2_proj_p_'
			self.Path = Paths.control_files+'test/L2_proj_p/'
		elif (self.name.find('L2_proj_h') != -1):
			self.name = 'Test_L2_proj_h_'
			self.Path = Paths.control_files+'test/L2_proj_h/'
		elif (self.name.find('linearization') != -1):
			self.name = 'Test_linearization_'
			self.Path = Paths.control_files+'test/linearization/'
		elif (self.name.find('Poisson') != -1):
			if (self.name.find('test') != -1):
				self.name = 'Test_Poisson_'
				self.Path = Paths.control_files+'test/Poisson/'
			else:
				self.name = 'Poisson'
				print("name:",self.name)
				EXIT_TRACEBACK()
		elif (self.name.find('Euler') != -1):
			if (self.name.find('test') != -1):
				self.name = 'Test_Euler_'
				self.Path = Paths.control_files+'test/Euler/'
			else:
				print("name:",self.name)
				EXIT_TRACEBACK()
		elif (self.name.find('NavierStokes') != -1):
			if (self.name.find('test') != -1):
				self.name = 'Test_NavierStokes_'
				self.Path = Paths.control_files+'test/NavierStokes/'
			else:
				print("name:",self.name)
				EXIT_TRACEBACK()
		else:
			print("name:",self.name)
			EXIT_TRACEBACK()

	def get_geo_dependencies(self):
		MeshTypes = self.MeshTypes
		for i in range(0,len(MeshTypes)):
			if (self.GeoDeps.find(MeshTypes[i].InputName) == -1):
				self.GeoDeps += MeshTypes[i].InputName + ' '

	def get_mesh_outputs(self):
		MeshTypes = self.MeshTypes
		for i in range(0,len(MeshTypes)):
			if (self.MeshOutputs.find(MeshTypes[i].InputName) == -1):
				self.MeshOutputs += MeshTypes[i].OutputName_from_meshesROOT + ' '


class MeshType_class:
	def __init__(self,name,prefix):
		self.name   = name
		self.prefix = prefix

		self.PDEName       = ''
		self.PDESpecifier  = ''
		self.Geometry      = ''
		self.GeomSpecifier = ''
		self.dim           = ''
		self.MeshLevel     = ''
		self.MeshCurving   = ''

		self.InputName  = ''
		self.OutputDir  = ''
		self.OutputName = ''
		self.OutputName_from_meshesROOT = ''

	def set_parameters(self,TestCase,Paths):
		fName = TestCase.Path+TestCase.name+self.prefix+self.name+'.ctrl'

		with open(fName) as f:
			for line in f:
				if ('PDEName' in line):
					self.PDEName = line.split()[1]
				if ('PDESpecifier' in line):
					self.PDESpecifier = line.split()[1]
				if ('Geometry' in line):
					self.Geometry = line.split()[1]
					self.InputName = line.split()[2]
				if ('GeomSpecifier' in line):
					self.GeomSpecifier = line.split()[1]
				if ('MeshCurving' in line):
					self.MeshCurving = line.split()[1]
				if ('Dimension' in line):
					self.dim = line.split()[1]
				if ('MeshLevel' in line):
					self.MeshLevel = line.split()[1]

		self.OutputDir  = self.Geometry + '/' + self.PDEName + '/'
		if (self.PDESpecifier.find("NONE") == -1):
			self.OutputDir += self.PDESpecifier + '/'
		if (self.GeomSpecifier.find("NONE") == -1):
			self.OutputDir += self.GeomSpecifier + '/'
		self.OutputName = self.OutputDir + self.Geometry + self.dim + 'D_'
		self.OutputDir  = Paths.meshes + self.OutputDir
		
		if (self.MeshCurving.find('Straight') != 0):
			self.OutputName += self.MeshCurving
		self.OutputName += self.name + self.MeshLevel + 'x.msh'

		self.OutputName_from_meshesROOT = self.OutputName
		self.OutputName = Paths.meshes + self.OutputName
