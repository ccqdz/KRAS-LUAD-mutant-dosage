# Release and DOI checklist

1. Create a public GitHub repository, suggested name:
   `KRAS-LUAD-mutant-dosage`
2. Upload the contents of this directory to the repository root.
3. Replace `REPOSITORY_URL_TO_BE_INSERTED` in `CITATION.cff` with the final GitHub URL.
4. Select and add the institution-approved software license.
5. Create a GitHub release/tag `v1.0.0` corresponding to the submitted manuscript.
6. Link the GitHub repository to Zenodo and enable release archiving.
7. Publish/re-publish GitHub release `v1.0.0`; Zenodo will mint a version DOI and a concept DOI.
8. Insert the concept DOI (preferred for the evolving repository) or the version DOI (preferred for the exact archived submission) into the manuscript Code availability statement.
9. Freeze the exact release used for submission; do not overwrite or force-update the archived tag.

Suggested manuscript wording after URLs are assigned:

> Code availability: The consolidated R workflow, minimum input manifest, software environment, and locked reproducibility audits are publicly available at [GITHUB_URL]. The archived release corresponding to this manuscript is available at Zenodo: https://doi.org/[ZENODO_DOI]. Patient-level controlled-access data are not redistributed and must be obtained from the originating resources under their applicable access conditions.
