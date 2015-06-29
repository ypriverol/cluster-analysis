__author__ = 'yperez'


class ClusterEntry:

    def __init__(self,
                 cluster_id,
                 precursor_mz,
                 size,
                 max_ratio,
                 max_il_ratio,
                 precursor_mz_range,
                 max_sequence,
                 max_sequence_count,
                 max_sequence_projects,
                 max_sequence_contaminant,
                 second_max_sequence,
                 second_max_sequence_count,
                 second_max_sequence_projects,
                 second_max_sequence_contaminant,
                 third_max_sequence,
                 third_max_sequence_count,
                 third_max_sequence_projects,
                 third_max_sequence_contaminant,
                 project_count,
                 assay_count,
                 species):

        self.clusterID                       = cluster_id
        self.precursor_mz                    = precursor_mz
        self.size                            = size
        self.max_ratio                       = max_ratio
        self.max_il_ratio                    = max_il_ratio
        self.precursor_mz_range              = precursor_mz_range
        self.max_sequence                    = max_sequence
        self.max_sequence_count              = max_sequence_count
        self.max_sequence_projects           = max_sequence_projects
        self.max_sequence_contaminant        = max_sequence_contaminant
        self.second_max_sequence             = second_max_sequence
        self.second_max_sequence_count       = second_max_sequence_count
        self.second_max_sequence_projects    = second_max_sequence_projects
        self.second_max_sequence_contaminant = second_max_sequence_contaminant
        self.third_max_sequence              = third_max_sequence
        self.third_max_sequence_count        = third_max_sequence_count
        self.third_max_sequence_projects     = third_max_sequence_projects
        self.third_max_sequence_contaminant  = third_max_sequence_contaminant
        self.project_count                   = project_count
        self.assay_count                     = assay_count
        self.species                         = species

    def clusterID (self) :
        return self.clusterID

    def precursorMz(self):
        return self.precursor_mz

    def size(self):
        return self.size

    def maxRatio(self):
        return self.max_ratio

    def maxIlRatio(self):
        return self.max_il_ratio

    def precursorMzRange(self):
        return self.precursor_mz_range

    def maxSequence(self):
        return self.max_sequence

    def maxSequenceCount(self):
        return self.max_sequence_count

    def maxSequenceProjects(self):
        return self.max_sequence_projects

    def maxSequenceContaminant(self):
        return self.max_sequence_contaminant

    def secondMaxSequence(self):
        return self.second_max_sequence

    def secondMaxSequenceCount(self):
        return self.second_max_sequence_count

    def secondMaxSequenceProjects(self):
        return self.second_max_sequence_projects

    def secondMaxSequenceContaminant(self):
        return self.second_max_sequence_contaminant

    def thirdMaxSequence(self):
        return self.third_max_sequence

    def thirdMaxSequenceCount(self):
        return self.third_max_sequence_count

    def thirdMaxSequenceProjects(self):
        return self.third_max_sequence_projects

    def thirdMaxSequenceContaminant(self):
        return self.third_max_sequence_contaminant

    def projectCount(self):
        return self.project_count

    def assayCount(self):
        return self.assay_count

    def species(self):
        return self.species