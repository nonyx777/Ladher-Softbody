extends Label


func _process(_delta: float):
	text = "FPS: " + str(Engine.get_frames_per_second())
