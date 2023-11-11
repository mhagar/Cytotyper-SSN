from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired
from wtforms import SubmitField, IntegerField


class TestForm(FlaskForm):
    user_upload = FileField(validators=[FileRequired()])
    user_upload_query = FileField(validators=[FileRequired()])
    user_min_occurences = IntegerField()
    user_min_gap_size = IntegerField()
    submit = SubmitField('GO!!')
