"""Detect execution identity without submitting a cloud job."""
import io
import unittest
from unittest.mock import Mock, patch
from urllib.error import URLError

class IdentityDefaultTests(unittest.TestCase):
    def test_explicit_values_do_not_contact_metadata(self):
        from scripts.glow.controller import resolve_execution_identity
        reader = Mock(side_effect=AssertionError("Metadata must not be called"))
        self.assertEqual(resolve_execution_identity("project-one", "worker@project-one.iam.gserviceaccount.com",
                                                   metadata_reader=reader),
                         ("project-one", "worker@project-one.iam.gserviceaccount.com"))
        reader.assert_not_called()

    def test_missing_values_use_the_attached_task_identity(self):
        from scripts.glow.controller import resolve_execution_identity
        values = {"project/project-id": "terra-project-one",
                  "instance/service-accounts/default/email": "pet-worker@terra-project-one.iam.gserviceaccount.com"}
        reader = Mock(side_effect=values.__getitem__)
        self.assertEqual(resolve_execution_identity(None, None, metadata_reader=reader),
                         ("terra-project-one", "pet-worker@terra-project-one.iam.gserviceaccount.com"))
        self.assertEqual([c.args[0] for c in reader.call_args_list], list(values))

    def test_each_override_can_be_set_separately(self):
        from scripts.glow.controller import resolve_execution_identity
        reader = Mock(return_value="pet-worker@terra-project-one.iam.gserviceaccount.com")
        project, account = resolve_execution_identity("other-project", None, metadata_reader=reader)
        self.assertEqual(project, "other-project")
        self.assertEqual(account, reader.return_value)
        reader.assert_called_once_with("instance/service-accounts/default/email")
        reader = Mock(return_value="terra-project-one")
        self.assertEqual(resolve_execution_identity(None, "worker@other-project.iam.gserviceaccount.com",
                                                    metadata_reader=reader),
                         ("terra-project-one", "worker@other-project.iam.gserviceaccount.com"))
        reader.assert_called_once_with("project/project-id")

    def test_detection_failure_requires_an_override_instead_of_another_account(self):
        from scripts.glow.controller import resolve_execution_identity
        for failure in (URLError("metadata unavailable"), ValueError("invalid metadata response")):
            with self.subTest(failure=failure):
                with self.assertRaisesRegex(ValueError, "--project"):
                    resolve_execution_identity(None, None, metadata_reader=Mock(side_effect=failure))
        with self.assertRaisesRegex(ValueError, "--service-account"):
            resolve_execution_identity("project-one", None, metadata_reader=Mock(return_value=""))

    def test_metadata_uses_google_header_short_timeout_and_no_proxy(self):
        from scripts.glow.controller import read_instance_metadata
        response = io.BytesIO(b"terra-project-one\n")
        response.headers = {"Metadata-Flavor": "Google"}
        opener = Mock()
        opener.open.return_value = response
        with patch("scripts.glow.controller.build_opener", return_value=opener) as factory:
            self.assertEqual(read_instance_metadata("project/project-id"), "terra-project-one")
        self.assertEqual(factory.call_args.args[0].proxies, {})
        request = opener.open.call_args.args[0]
        self.assertEqual(request.full_url, "http://metadata.google.internal/computeMetadata/v1/project/project-id")
        self.assertEqual(request.get_header("Metadata-flavor"), "Google")
        self.assertEqual(opener.open.call_args.kwargs["timeout"], 3)

    def test_metadata_response_must_come_from_google(self):
        from scripts.glow.controller import read_instance_metadata
        response = io.BytesIO(b"project-one")
        response.headers = {}
        opener = Mock()
        opener.open.return_value = response
        with patch("scripts.glow.controller.build_opener", return_value=opener):
            with self.assertRaisesRegex(ValueError, "metadata"):
                read_instance_metadata("project/project-id")
