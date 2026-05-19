#' RDoC Term-to-Domain Reference
#'
#' Returns the internal RDoC term list used by [rdoc_decode()].
#'
#' @return A data frame with columns `Domain` and `Term`.
#' @export
rdoc_terms_reference <- function() {
  data.frame(
    Domain = c(
      "AR", "AR", "CS", "CS", "CS", "CS", "CS", "CS", "CS",
      "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS", "CS",
      "CS", "NV", "NV", "NV", "PV", "PV", "PV", "PV", "PV",
      "SP", "SP", "SP", "SP", "SS", "SS", "SS", "SS", "SS", "SS"
    ),
    Term = c(
      "Arousal", "Circadian Rhythm", "Attention",
      "CC - Goal Selection (Goal)", "CC - Goal Selection (Selection)",
      "CC - Performance Monitoring", "CC - Response Selection",
      "CC - Inhibition/Suppression", "Declarative Memory",
      "Declarative Memory/WM", "Language", "Perception - Auditory",
      "Perception - Multimodal", "Perception - Olfactory", "Perception - Somatosensory",
      "Perception - Visual", "WM - Interference Control", "WM - Active Maintenance",
      "WM - Flexible Updating", "Acute Threat ('Fear')", "Loss", "Potential Threat ('Anxiety')",
      "RL - Reward Prediction Error", "RL - Probabilistic and Reinforcement Learning",
      "Reward Anticipation", "Reward Responsiveness", "Reward Valuation",
      "Affiliation and Attachment", "Animacy Perception", "Action Perception",
      "Social Communication", "Agency and Ownership", "Motor Execution",
      "Motor Actions", "Motor Inhibition and Termination", "Motor Initiation",
      "Innate Motor Patterns"
    ),
    term_id = c(
      "arousal", "circadian_sleep_wakefulness", "attention",
      "goal", "selection",
      "performance_monitoring", "response_selection",
      "suppression", "declarative_memory",
      "memory", "language", "perception_auditory",
      "perception_multimodal", "perception_olfactory", "somatosensory_perception",
      "visual_perception", "interference_control",
      "working_memory_maintenance_working_memory_capacity",
      "working_memory_updating", "fear", "loss", "anxiety",
      "prediction_error", "probabilistic_learning_reinforcement_learning",
      "reward_anticipation", "reward_response", "reward_probability", 
      "attachment", "animacy", "action_perception",
      "communication", "sensory_agency_sensory_ownership", "execution",
      "motor_action", "motor_inhibition", "motor_initiation",
      "motor_pattern"
    ),
    stringsAsFactors = FALSE
  )
}

#' Read Example RDoC Correlation Data
#'
#' Loads the bundled TSV example data with columns `Domain`, `Term`, `r`, and `p`.
#'
#' @return A data frame.
#' @export
rdoc_example_data <- function() {
  path <- system.file(
    "extdata",
    "methylphenidate_mean_RDoC_absolute-values.tsv",
    package = "rdocodeR"
  )
  if (identical(path, "")) {
    stop("Bundled example file is missing from this installation.", call. = FALSE)
  }
  utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

#' Path to Bundled RDoC Term Maps
#'
#' Returns the path of the bundled `rdoc_terms.fsavg6.lh.rh.rds` file.
#'
#' @return Absolute file path to the bundled term map RDS file.
#' @export
rdoc_terms_file <- function() {
  path <- system.file("extdata", "rdoc_terms.fsavg6.lh.rh.rds", package = "rdocodeR")
  if (identical(path, "")) {
    stop("Bundled RDoC term map file is missing from this installation.", call. = FALSE)
  }
  path
}
