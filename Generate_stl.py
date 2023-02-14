# -*- coding:UTF-8 -*-
import rhinoscriptsyntax as rs
import math
import os
import shutil
Mesh = rs.GetObject("请选择模型")
filename = rs.DocumentName()
filepath = rs.DocumentPath()
source_system='F:\\课题研究\\抗风设计\\案例\\绍兴植物园\\files\\system\\'
source_constant = 'F:\\课题研究\\抗风设计\\案例\\绍兴植物园\\files\\constant\\'
source_0 = 'F:\\课题研究\\抗风设计\\案例\\绍兴植物园\\files\\0'
wind_pressure = rs.GetReal("请输入基本风压(kN/m2)")
Site_Type = rs.GetString("请输入粗糙度类别(A,B,C,D)")
#wind_pressure = 0.5
#Site_Type = "A"
degrees = 0    #每次转动的角度增量，单位度数
Num = 1        #生成模型数量
basic_velocity = math.sqrt(wind_pressure*1600)
if Mesh is not None and rs.IsMesh(Mesh):
    points = rs.MeshVertices(Mesh)
    xMin = 0
    xMax = 0
    yMin = 0
    yMax = 0
    zMin = 0
    zMax = 0
    for point in points:
        if point[0]<xMin:
            xMin = point[0]
        if point[0]>xMax:
            xMax = point[0]     
        if point[1]<yMin:
            yMin = point[1]
        if point[1]>yMax:
            yMax = point[1]                   
        if point[2]<zMin:
            zMin = point[2]
        if point[2]>zMax:
            zMax = point[2]
    CP = rs.MeshAreaCentroid(Mesh)
    translation = -CP
    rs.MoveObject(Mesh,(translation[0],translation[1],-zMin))

    for j in range(0,Num):
        radians = math.radians(degrees)
        c = math.cos(radians)
        s = math.sin(radians)
        matrix = []
        matrix.append( [c,-s, 0, 0] )
        matrix.append( [s, c, 0, 0] )
        matrix.append( [0, 0, 1, 0] )
        matrix.append( [0, 0, 0, 1] )
        rs.TransformObject( Mesh, matrix )
        model_filename = "WindModel_"+str(int(degrees*j))
        if not os.path.exists(model_filename+'\\constant'):
            os.mkdir(model_filename+'\\constant')
            if not os.path.exists(model_filename+'\\constant\\triSurface'):
                os.mkdir(model_filename+'\\constant\\triSurface')
        shutil.copy(source_constant+'transportProperties', model_filename+'\\constant')
        shutil.copy(source_constant+'turbulenceProperties', model_filename+'\\constant')
        SavePath = filepath+model_filename+"\\constant\\triSurface\\" +"WindModel_"+str(int(degrees*j))+".stl"
        rs.UnselectAllObjects()
        rs.SelectObject(Mesh)
        comm_opts='_Version=7 _SaveSmall=_No _GeometryOnly=_No _SaveTextures=_No _SavePlugInData=_No'
        rs.Command("_-Export " +comm_opts + " "+ SavePath+ " _Enter" +" _Enter")
        if not os.path.exists(model_filename+'\\system'):
            os.mkdir(model_filename+'\\system')
    #Write blocMeshDict
        Box1 = max(xMax-xMin,yMax-yMin)
        Box2 = zMax-zMin
        Box_xMin = (xMin-2*Box1)
        Box_xMax = (xMax+2*Box1)
        Box_yMin = (yMin-1*Box1)
        Box_yMax = (yMax+1*Box1)
        Box_zMin = 0
        Box_zMax = 2*max(Box1,Box2)
        with open(source_system+"blockMeshDict",'r') as block_file:
            lines = block_file.readlines()
            block_file.close
            lines[24] = "  ( {} {} {} )\n".format(Box_xMin,Box_yMin,0)
            lines[25] = "  ( {} {} {} )\n".format(Box_xMax,Box_yMin,0)
            lines[26] = "  ( {} {} {} )\n".format(Box_xMax,Box_yMax,0)
            lines[27] = "  ( {} {} {} )\n".format(Box_xMin,Box_yMax,0)
            lines[28] = "  ( {} {} {} )\n".format(Box_xMin,Box_yMin,Box_zMax)
            lines[29] = "  ( {} {} {} )\n".format(Box_xMax,Box_yMin,Box_zMax)
            lines[30] = "  ( {} {} {} )\n".format(Box_xMax,Box_yMax,Box_zMax)
            lines[31] = "  ( {} {} {} )\n".format(Box_xMin,Box_yMax,Box_zMax)
        with open(model_filename+"\\system\\blockMeshDict",'w+') as block_file:
            for line in lines:
                block_file.write(line)
            block_file.close
    #Write snappyHexMeshDict
        with open(source_system+"snappyHexMeshDict",'r') as snappyHexMesh_file:
            lines = snappyHexMesh_file.readlines()
            snappyHexMesh_file.close
            lines[32] = "    {}         {{ type triSurfaceMesh; name Model_wall; }}\n".format("WindModel_"+str(int(degrees*j))+".stl")
            lines[33] = "    MeshRefinementBox_000 {{ type searchableBox; min ({} {} {}); max ({} {} {}); }}\n".format(0.5*Box_xMin,0.5*Box_yMin,0.5*Box_zMin,0.5*Box_xMax,0.5*Box_yMax,0.5*Box_zMax)
            lines[34] = "    MeshRefinementBox_001 {{ type searchableBox; min ({} {} {}); max ({} {} {}); }}\n".format(0.3*Box_xMin,0.3*Box_yMin,0.3*Box_zMin,0.3*Box_xMax,0.3*Box_yMax,0.3*Box_zMax)
            lines[101] = "       Model_wall  { level (2 4); patchInfo { type wall; }}\n"
            lines[142] = "    locationInMesh ({} {} {});\n".format(0.5*(Box_xMin+Box_xMax),0.5*(Box_yMin+Box_yMax),0.5*(Box_zMin+Box_zMax))
            lines[127] = "        Model_wall\n"
            lines[233] = "        Model_wall\n"
        with open(model_filename+"\\system\\snappyHexMeshDict",'w+') as snappyHexMesh_file:
            for line in lines:
                snappyHexMesh_file.write(line)
            snappyHexMesh_file.close
    #Write decomposeParDict
        with open(source_system+"decomposeParDict",'r') as decomposePar_file:
            lines = decomposePar_file.readlines()
        with open(model_filename+"\\system\\decomposeParDict",'w+') as decomposePar_file:
            for line in lines:
                decomposePar_file.write(line)
            decomposePar_file.close
    #Write controlDict
        with open(source_system+"controlDict",'r') as control_file:
            lines = control_file.readlines()
        with open(model_filename+"\\system\\controlDict",'w+') as control_file:
            for line in lines:
                control_file.write(line)
            control_file.close
    #Write fvSolution
        with open(source_system+"fvSolution",'r') as fvSolution_file:
            lines = fvSolution_file.readlines()
        with open(model_filename+"\\system\\fvSolution",'w+') as fvSolution_file:
            for line in lines:
                fvSolution_file.write(line)
            fvSolution_file.close
    #Write fvSchemes
        with open(source_system+"fvSchemes",'r') as fvSchemes_file:
            lines = fvSchemes_file.readlines()
        with open(model_filename+"\\system\\fvSchemes",'w+') as fvSchemes_file:
            for line in lines:
                fvSchemes_file.write(line)
            fvSchemes_file.close
    #Write changeDictionary
        with open(source_system+"changeDictionaryDict",'r') as changeDictionary_file:
            lines = changeDictionary_file.readlines()
        with open(model_filename+"\\system\\changeDictionaryDict",'w+') as changeDictionary_file:
            for line in lines:
                changeDictionary_file.write(line)
            changeDictionary_file.close
    #Write initial Condition           
        if not os.path.exists(model_filename+"\\0"):
            os.mkdir(model_filename+"\\0")
        for root, dirs, files in os.walk(source_0):
            for file in files:
                src_file = os.path.join(root, file)
                shutil.copy(src_file, model_filename+"\\0")
        I = 0.01 #湍流强度默认取1%               
        k = 1.5*(basic_velocity*I)**2
        L = 0.1643*max(yMax-yMin,zMax-zMin)
        epsilon = 0.09*k**1.5/L    #0.09为默认常数
        
        with open(model_filename+"\\0\\k",'r+') as k_file:
            lines = k_file.readlines()
            lines[23] = "internalField uniform {};\n".format(k)
            lines[30] = "        value uniform {};\n".format(k)
            k_file.seek(0,0)
            for line in lines:
                k_file.write(line)
            k_file.close               

        with open(model_filename+"\\0\\epsilon",'r+') as epsilon_file:
            lines = epsilon_file.readlines()
            lines[20] = "internalField uniform {};\n".format(epsilon)
            lines[27] = "        value uniform {};\n".format(epsilon)
            epsilon_file.seek(0,0)
            for line in lines:
                epsilon_file.write(line)
            epsilon_file.close    


    #Caculate wind profile
        Height_z = [5*x for x in range(0, 21, 1)]
        Height_z += [100+50*x for x in range(1, 11, 1)]
        parameter = []
        parameter.append( [0.12,0.15, 0.22, 0.3] )
        parameter.append( [1.133,1, 0.738, 0.512] )
        parameter.append( [5,10, 15, 20] )
        parameter.append( [300,350, 450, 550] )
        if Site_Type == "A":
            rough_facter = parameter[0][0]
            k =  parameter[1][0]
            H1 = parameter[2][0]
            H2 = parameter[3][0]
        elif Site_Type == "B":
            rough_facter = parameter[0][1]
            k =  parameter[1][1]
            H1 = parameter[2][1]
            H2 = parameter[3][1]
        elif Site_Type == "C":
            rough_facter = parameter[0][2]
            k =  parameter[1][2]
            H1 = parameter[2][2]
            H2 = parameter[3][2]
        elif Site_Type == "D":
            rough_facter = parameter[0][3]
            k =  parameter[1][3]
            H1 = parameter[2][3]
            H2 = parameter[3][3]
        U = [ k*basic_velocity*(min(max(H1,z),H2)/10)**(rough_facter) for z in Height_z]
        U = [ round(x,3) for x in U]
        #Write wind profile
        with open(filepath+model_filename+"\\UProfile.txt",'w+') as UProfile_file:
            UProfile_file.write("(\n")
            for i,z in enumerate(Height_z):
                UProfile_file.write("({} ({} 0 0))\n".format(z,U[i]))
            UProfile_file.write(")")
            UProfile_file.close