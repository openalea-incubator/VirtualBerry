"""This is an example on running HydroShoot on a potted grapevine with a
simple shoot architecture.
"""

from json import load, dump
from pathlib import Path

from openalea.mtg import traversal, mtg
from openalea.mtg.mtg import MTG
from openalea.plantgl.all import Scene
from hydroshoot import architecture, display
from hydroshoot import io, model, initialisation


def build_mtg(path_file: Path, is_show_scene: bool = True) -> (mtg.MTG, Scene):
    grapevine_mtg = architecture.vine_mtg(file_path=path_file)

    for v in traversal.iter_mtg2(grapevine_mtg, grapevine_mtg.root):
        architecture.vine_phyto_modular(grapevine_mtg, v)
        architecture.vine_mtg_properties(grapevine_mtg, v)
        architecture.vine_mtg_geometry(grapevine_mtg, v)
        architecture.vine_transform(grapevine_mtg, v)

    # Display of the plant mock-up (result in 'fig_01_plant_mock_up.png')
    mtg_scene = display.visu(grapevine_mtg, def_elmnt_color_dict=True, scene=Scene(), view_result=is_show_scene)
    return grapevine_mtg, mtg_scene


def preprocess_inputs(grapevine_mtg: MTG, path_project_dir: Path, psi_soil: float, gdd_since_budbreak: float,
                      scene: Scene):
    inputs = io.HydroShootInputs(g=grapevine_mtg, path_project=path_project_dir, scene=scene, psi_soil=psi_soil,
                                 gdd_since_budbreak=gdd_since_budbreak)
    io.verify_inputs(g=grapevine_mtg, inputs=inputs)
    grapevine_mtg = initialisation.init_model(g=grapevine_mtg, inputs=inputs)

    static_data = {'form_factors': {s: grapevine_mtg.property(s) for s in ('ff_sky', 'ff_leaves', 'ff_soil')}}
    static_data.update({'Na': grapevine_mtg.property('Na')})

    dynamic_data = {}
    inputs_hourly = io.HydroShootHourlyInputs(psi_soil=inputs.psi_soil_forced, sun2scene=inputs.sun2scene)
    for date_sim in inputs.params.simulation.date_range:
        inputs_hourly.update(
            g=grapevine_mtg, date_sim=date_sim, hourly_weather=inputs.weather[inputs.weather.index == date_sim],
            psi_pd=inputs.psi_pd, params=inputs.params)

        grapevine_mtg, diffuse_to_total_irradiance_ratio = initialisation.init_hourly(
            g=grapevine_mtg, inputs_hourly=inputs_hourly, leaf_ppfd=inputs.leaf_ppfd,
            params=inputs.params)

        dynamic_data.update({grapevine_mtg.date: {
            'diffuse_to_total_irradiance_ratio': diffuse_to_total_irradiance_ratio,
            'Ei': grapevine_mtg.property('Ei'),
            'Eabs': grapevine_mtg.property('Eabs')}})

    path_preprocessed_inputs = path_project_dir / 'preprocessed_inputs'
    with open(path_preprocessed_inputs / 'static.json', mode='w') as f_prop:
        dump(static_data, f_prop)
    pass
    with open(path_preprocessed_inputs / 'dynamic.json', mode='w') as f_prop:
        dump(dynamic_data, f_prop)
    pass


if __name__ == '__main__':
    path_project = Path(__file__).parent
    path_preprocessed_data = path_project / 'preprocessed_inputs'

    g, scene = build_mtg(path_file=path_project / 'grapevine_pot.csv', is_show_scene=False)
    # preprocess_inputs(grapevine_mtg=g, path_project_dir=path_project, psi_soil=-0.5, gdd_since_budbreak=1000., scene=scene)

    with open(path_preprocessed_data / 'static.json') as f:
        static_inputs = load(f)
    with open(path_preprocessed_data / 'dynamic.json') as f:
        dynamic_inputs = load(f)

    summary_results = model.run(g=g, wd=path_project, scene=scene, gdd_since_budbreak=1000.,
                                form_factors=static_inputs['form_factors'],
                                leaf_nitrogen=static_inputs['Na'],
                                leaf_ppfd=dynamic_inputs,
                                path_output=path_project / 'output' / 'time_series_with_preprocessed_data.csv')
