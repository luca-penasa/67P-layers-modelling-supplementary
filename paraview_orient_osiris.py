### Paraview macro to align the camera with any OSIRIS image

# numpy for arrays manipulations
import numpy as np
# we use spiceypy, other toolkits might be used
import spiceypy
import os, sys
import re
import vtk

# the location where Rosetta spice kernels are kept.
# the kernels are public so you can find them here:
# ftp://anonymous@naif.jpl.nasa.gov/pub/naif/ROSETTA
# an example:
spice_kernels_dir = "/home/luca/spice_kernels"


# you also need a furnsh file. Have a look at the SPICE kernels documentation.
# its a list of the spice kernels needed for the computation
furnsh_file = "kernels_CG.furnsh"
furnsh_file = "ROS_V111.TM"

# the picture to align to
# modify the code accordingly if you just have the time string
fname = "NAC_2016-01-27T16.27.58.970Z_ID30_1397549500_F22"

# if screenshot is needed set screenshotfile to some path.
# e.g. screenshotfile = "/home/luca/test.png"
screenshotfile = None

def vtkmatrix4x4_to_array(vtkmat):
    '''
    transform a vtk matrix to a numpy array
    :param vtkmat: the vtkMatrix4x4 object
    :return: the array
    '''
    scipy_array = np.zeros((4, 4), 'd')
    # we get the elements one at the time
    for i in range(0, 4):
        for j in range(0, 4):
            scipy_array[i][j] = vtkmat.GetElement(i, j)
    return scipy_array

def extractTimeString(name):
    '''
    From a filename to the time string by using regex
    Notice that time encoded in OSIRIS pictures is not granted to match the exact time the picture was taken.
    For this it would be sensitive to use the time encoded in the pictures metadata or convert the spacecraft clock to
    a date/time.
    :param name: an OSIRIS image name e.g.  def getNacCenterAndRotationAtTime(imname):"NAC_2014-08-06T02.43.16.574Z_ID30."
    :return: only the time substring 2014-08-06T02.43.16.574
    '''
    m = re.search("[A-Z]{3}_(.*)T(.*)[A-Z]{1}_.*", name)
    parta = m.group(1);
    partb = m.group(2);
    partb = partb.replace(".", ":")
    partb = list(partb)
    partb[-4] = '.'
    partb = ''.join(partb)
    return parta + "T" + partb;


def getNacCenterAndRotationAtTime(imname):
    '''
    Extract the rotation of the spacecraft and its position for the given picture
    Notice that WAV and NAC can be considered aligned on the same axis and centred in the spacecraft
    so it can be used the same way independently by the camera
    :param imname: the input picture filename
    :return: the centre and the rotation matrix of the spacecraft
    '''
    os.chdir(spice_kernels_dir) # or spyceypy will not work
    spiceypy.furnsh(spice_kernels_dir + os.path.sep + furnsh_file) #load the kernels from furnsh

    tstr = extractTimeString(imname) # get only the time string from the filename
    t = spiceypy.str2et(tstr) # transform the string to a time

    # get the position of Rosetta spacecraft in 67p refence frame, centred in 67P + light time corr (not really needed)
    pos = spiceypy.spkpos("ROSETTA", t, "67P/C-G_CK", "LT+S", "67P/C-G")[0]
    # pos = spiceypy.spkpos("ROSETTA", t, "ROS_OSIRIS_NAC", "LT+S", "67P/C-G")[0]

    print ("Found position: {}".format(pos))

    #and the rotation - we use NAC but it is the same for the WAC
    # R = spiceypy.pxform("67P/C-G_CK", "ROS_OSIRIS_NAC", t);
    # rotation matrix of the nac camera, in respect to 67P reference frame
    R = spiceypy.pxform("ROS_OSIRIS_NAC", "67P/C-G_CK", t);
    return np.array(pos), np.array(R)

def whatCameraIs(filename):
    '''
    Match if the camera is NAC or WAC from the filename
    :param filename: the input filename of the picture
    :return: WAC, NAC or None
    '''
    if re.match(".*NAC.*", filename):
        return "NAC"
    elif re.match(".*WAC.*", filename):
        return "WAC"
    else:
        return None


if whatCameraIs(fname) is not None:
    # and the active camera
    camera = GetActiveCamera()
    # ensure parallel projection is off
    view = GetActiveView()
    view.CameraParallelProjection = 0

    # get centre and rotation of camera
    c, R = getNacCenterAndRotationAtTime(fname)

    # we create a homogeneous transformation matrix,
    # embedding centre and rotation
    Rt = np.zeros((4,4))
    Rt[3,3] = 1
    Rt[:3,:3] = R
    Rt[:3,3] = c

    #we do a total reset to start in a known position
    camera.SetFocalPoint(0, 0, 1)
    camera.SetPosition(0, 0, 0)
    camera.SetViewUp((0, 1, 0))
    camera.SetRoll(-90)
    camera.SetWindowCenter(0, 0)

    # create a transform
    t = vtk.vtkTransform()

    #finally set to the transform
    t.SetMatrix(Rt.ravel().tolist())

    # and apply the transform the camera
    camera.ApplyTransform(t)

    # set the FoV depending on the camera
    if whatCameraIs(fname) == "NAC":
        camera.SetViewAngle(2.35)
    else:
        camera.SetViewAngle(12)  # 12 for WAC

    Render()

    # we might save a screenshot
    if screenshotfile is not None:
        SaveScreenshot(screenshotfile, ImageQuality=100, viewOrLayout=view, TransparentBackground=True,
                   ImageResolution=(2048, 2048))

