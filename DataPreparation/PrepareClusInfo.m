function Params = PrepareClusInfo(KiloSortPaths, Params, RawDataPathsInput)

if nargin<3
Params = ExtractKilosortData(KiloSortPaths, Params)

else
Params = ExtractKilosortData(KiloSortPaths, Params, RawDataPathsInput)
end