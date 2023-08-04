import argparse, json, os

default_dict = {
    "opt_numConfs": 5,
    "opt_numThreads": 0,
    "opt_maxIters": 1500,
    "drawing_subImgSize_edge": 250,
    "drawing_default": "polymer.png",
    "MV_gridSpacing": 0.2,
    "MV_boxMargin": 2.0,
    "plot_dataPoint": "o",
    "plot_Filename": "Size-dependent-stats.png",
    "dielectricModel": 2,
    "dielectricConstant": 78,
    "NB_THRESH": 100
}


def readJson(filepath):
    with open(filepath, "r") as S: #this is where many defaults are set so they can easily be changed.
        settings_dict = json.load(S)
    return settings_dict


def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--show', default = False, action = "store_true", help = "prints all settings")
    parser.add_argument('-w', '--write', default = False, action = "store_true", help = "writes a copy of settings.json to cwd for user to edit manually.")
    args = parser.parse_args()
    return args


def writeJson(dict, path):
    json_object = json.dumps(dict, indent=4)
    with open(path, "w") as outfile:
        outfile.write(json_object)


def main():
    args = getArgs()
    settingsFile = "HXSettings.json"
    if args.show:
        if os.path.exists(settingsFile):
            print(f"All hydrophobicity_explorer packages will use these user-defined settings in {settingsFile}")
            user_settings = readJson(settingsFile)
            print(user_settings)
        else:
            print(default_dict)
    if args.write:
        writeJson(default_dict, settingsFile)
        print(f"Wrote {settingsFile}")


if __name__ == "__main__":
    main()
