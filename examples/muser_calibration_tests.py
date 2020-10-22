import logging
import sys

import numpy
from matplotlib import pylab

pylab.rcParams['figure.figsize'] = (8.0, 8.0)
pylab.rcParams['image.cmap'] = 'rainbow'

import matplotlib.pyplot as plt

import astropy.constants as constants

from rascil.data_models import PolarisationFrame
from rascil.processing_components import create_blockvisibility_from_ms, \
    advise_wide_field, create_image_from_visibility, convert_blockvisibility_to_stokesI, \
    gaintable_plot

from rascil.processing_components.flagging.operations import flagging_blockvisibility, \
    flagging_blockvisibility_with_bl

from rascil.workflows.rsexecute.execution_support import rsexecute, get_dask_client
from rascil.workflows import predict_list_rsexecute_workflow
from rascil.processing_components import create_calibration_controls
from rascil.workflows import calibrate_list_rsexecute_workflow

log = logging.getLogger("logger")
log.setLevel(logging.DEBUG)
log.addHandler(logging.StreamHandler(sys.stdout))

mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

if __name__ == "__main__":
    
    imaging_context = 'ng'

    client = get_dask_client(create_cluster=False, n_workers=4,
                             threads_per_worker=4, memory_limit="128GB")
    rsexecute.set_client(use_dask=False, client=client)
    
    # %%
    
    channels = range(60 - 4, 60 + 4)
    
    bvis_list = [rsexecute.execute(create_blockvisibility_from_ms)
                 ("/Users/wangfeng/work/muserdata/output/selfmodel.ms", start_chan=chan, end_chan=chan)[0]
                 for chan in channels]
    #
    # Remove this to check out circularnp
    # bvis_list = [rsexecute.execute(convert_blockvisibility_to_stokesI)(bv)
    #              for bv in bvis_list]
    # image_stokes = PolarisationFrame("stokesI")
    
    image_stokes = PolarisationFrame("stokesIV")

    def flag(bvis):
        bad = numpy.abs(bvis.vis) > 6e4
        print(numpy.sum(bad))
        bvis.flags[...][bad] = 1
        
        bvis = flagging_blockvisibility(bvis, antenna=[8, 9, 10, 11, 27])
        do_baselines = False
        if do_baselines:
            baseline = [[4, 0], [4, 1], [5, 4], [21, 4], [24, 4], [25, 4], [26, 4], [27, 4], [28, 4], [29, 4], [30, 4],
                        [31, 4],
                        [32, 4], [36, 4], [38, 4], [39, 4]]
            baseline.append(
                [[17, 4], [17, 13], [19, 17], [26, 17], [27, 17], [28, 17], [29, 17], [30, 17], [31, 17], [39, 17]])
            baseline.append([[19, 0], [19, 2], [19, 3], [19, 5], [19, 6], [19, 8], [19, 9], [19, 10], [19, 11], [19, 17]])
            baseline.append(
                [[20, 19], [24, 19], [25, 19], [26, 19], [27, 19], [28, 19], [29, 19], [30, 19], [31, 19], [34, 19],
                 [36, 19], [37, 19], [39, 19]])
            for i in range(28):
                baseline.append([28, i])
            baseline.append([[30, 28], [32, 28], [33, 28], [34, 28], [35, 28], [36, 28], [37, 28], [38, 28], [39, 28]])
            for i in range(28):
                baseline.append([29, i])
            baseline.append([[32, 29], [33, 29], [34, 29], [35, 29], [36, 29], [37, 29], [38, 29], [39, 29]])
            for i in range(29):
                baseline.append([30, i])
            baseline.append([[32, 30], [33, 30], [34, 30], [35, 30], [36, 30], [37, 30], [38, 30], [39, 30]])
            bvis = flagging_blockvisibility_with_bl(bvis, baseline)
        return bvis
    
    
    # %%
    bvis_list = [rsexecute.execute(flag)(bv) for bv in bvis_list]
    bvis_list = rsexecute.compute(bvis_list, sync=True)

    print(bvis_list[0])
    advice = advise_wide_field(bvis_list[-1], guard_band_image=0.5)
    
    cellsize = advice["cellsize"]
    npixel = advice["npixels2"]

    def read_sun_disk_data():
        sundisk_file_path = "quietsundisk400_2000MHz.txt"
        sundisk_file = open(sundisk_file_path, "r")
        sundisk_dict = {}
        try:
            while True:
                linekey = sundisk_file.readline()
                linevalue = sundisk_file.readline()
                if linekey and linevalue:
                    linekey = str(int(float(linekey.strip())))
                    value = list(map(eval, linevalue.strip().split(',')[1:]))
                    npvalue = []
                    for i in range(200):
                        npvalue.append(numpy.float32(value[i]))
                    sundisk_dict[linekey] = numpy.array(npvalue)
                else:
                    return sundisk_dict````
        finally:
            sundisk_file.close()
    
    sundisk_dict = read_sun_disk_data()
    
    
    def fill_solar_model_jy(m):
        # fov - minute
        size = m.data.shape[2]
        cellsize = m.wcs.wcs.cdelt[1] * numpy.pi / 180.0
        fov = cellsize * 180 / numpy.pi * size * 60
        # sun disk = 32 arc minute,
        sun_radius = 16
        scale = sun_radius * 2 / fov
        for chan in range(m.data.shape[0]):
            frequency = str(int(round(m.frequency[chan] * 1e-6)))
            assert frequency in sundisk_dict.keys(), "Key {} not present".format(frequency)
            sun_disk = sundisk_dict[frequency]
            for i in range(size):
                for j in range(size):
                    radius = numpy.sqrt((i - size // 2) ** 2 + (j - size // 2) ** 2) * fov / size
                    if radius <= sun_radius:
                        # print(int(round(radius/32*200,0)))
                        m.data[chan, 0, i, j] = sun_disk[int(round(radius / 32 * 200))]
            # The second polarisation is Q which should be zero (or much less than I)
            # m.data[chan,1,...] = 0.0
            # Convert from brightness temperature to Jy
            wavelength = constants.c.value / m.frequency[chan]
            t_to_jy = 1e26 * 2 * constants.k_B.value * cellsize ** 2 / wavelength ** 2
            m.data *= t_to_jy
            return m

    model_list = [rsexecute.execute(create_image_from_visibility)
                  (bvis, cellsize=cellsize, npixel=npixel,
                   polarisation_frame=image_stokes)
                  for bvis in bvis_list]

    model_list = [rsexecute.execute(fill_solar_model_jy, nout=1)(m) for m in model_list]
    
    bvis_list = rsexecute.scatter(bvis_list)
    from rascil.workflows import zero_list_rsexecute_workflow
    
    model_bvis_list = zero_list_rsexecute_workflow(bvis_list)
    model_bvis_list = predict_list_rsexecute_workflow(model_bvis_list, model_list, context=imaging_context)
    model_bvis_list = rsexecute.compute(model_bvis_list, sync=True)
    model_bvis_list = rsexecute.scatter(model_bvis_list)
    
    controls = create_calibration_controls()
    controls['G']['first_selfcal'] = 0
    controls['G']['phase_only'] = False
    for t_sol in [0.01, 0.1, 1.0, 10.0, 200.0]:
        controls['G']['timeslice'] = t_sol
        cal_graph = calibrate_list_rsexecute_workflow(bvis_list,
                                                      model_bvis_list,
                                                      gt_list=None,
                                                      calibration_context='G',
                                                      controls=controls,
                                                      global_solution=False)
        cal_bvis_list, gt_list = rsexecute.compute(cal_graph, sync=True)
        
        plt.clf()
        for igt, gt in enumerate(gt_list):
            title="G_{:.2f}s_{:.3f}GHz".format(t_sol, cal_bvis_list[igt].frequency[0]/1e9)
            gaintable_plot(gt['G'], 'G', title=title)
            plt.savefig("gain_{}.png".format(title))
            plt.show(block=False)
