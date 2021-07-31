from kivy.app import App
from kivy.uix.label import Label

from kivy.core.text import LabelBase, DEFAULT_FONT
from kivy.resources import resource_add_path

resource_add_path('C:\Windows\Fonts')
LabelBase.register(DEFAULT_FONT, 'meiryo.ttc') #日本語が使用できるように日本語フォントを指定する

class TextApp(App):
    def build(self):
        return Label(text='みなさん、こんにちは、ゆっひです。')

TextApp().run()        